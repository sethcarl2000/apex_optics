//geometry headers
#include <TGeoManager.h> 
#include <TGeoMaterial.h> 
#include <TGeoMedium.h> 
#include <TGeoMatrix.h> 
#include <TGeoArb8.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TEveManager.h> 
#include <TEvePointSet.h> 
#include <TEveArrow.h>
#include <TBrowser.h>  
#include <SieveHole.h> 
#include <ApexOptics.h>
#include <vector>  
#include <map> 
#include <cmath> 
#include <TMath.h> 
#include <TVector3.h> 
#include <TCanvas.h> 

//TGMainFrame windows
#include <TGWindow.h> 
#include <TGFrame.h> 
#include <TRootEmbeddedCanvas.h>
#include <TGButton.h>
#include <TGClient.h>  

using namespace std; 
//this is an enum to track which objects are drawn 
namespace Objects {
    enum NameBit : char16_t {
        kNone       = 0, 
        kSieve_LHRS = 1 << 0,           //left sieve 
        kSieve_RHRS = 1 << 1,           //right sieve
        kBeam       = 1 << 2,           //beam 
        kTarget_production = 1 << 3,    //production foils
        kTrack_L_real      = 1 << 4,    // LHRS - actual track
        kTrack_L_fg        = 1 << 5,    // LHRS - first-guess reconstruction
        kTrack_L_spread    = 1 << 6     // LHRS - spread of possible track trajectories
    };
}
//the definition of our geometry class
//_______________________________________________________________________________________________________________________________
class GeometryFrame : public TGMainFrame {
private: 
    TGVerticalFrame *fCanvFrame; 
    TRootEmbeddedCanvas *fECanvas; 

    TGeoManager *fGeom; 

    Objects::NameBit fDrawnObjects = Objects::kNone; 

    //a map of buttons, whith their associated name-bit. 
    std::vector<std::map<Objects::NameBit, TGCheckButton*>> fButtons{}; 

    void CreateGeometry(); 

public: 
    GeometryFrame(const TGWindow* p, UInt_t w, UInt_t h, const bool is_RHRS=false); 
    ~GeometryFrame(); 

    ClassDef(GeometryFrame,1); 
};
//_______________________________________________________________________________________________________________________________

GeometryFrame::~GeometryFrame() 
{
    fGeom->~TGeoManager(); 
    Cleanup(); 
}

void GeometryFrame::CreateGeometry() 
{
   //noop  
}

//_______________________________________________________________________________________________________________________________
GeometryFrame::GeometryFrame(const TGWindow* p, UInt_t w, UInt_t h, const bool is_RHRS)
    : TGMainFrame(p, w, h)
{
    fGeom = new TGeoManager("simple1", "Simple geometry");
    
    //--- define some materials
    TGeoMaterial *mat_vacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *mat_Al     = new TGeoMaterial("Al", 26.98, 13, 2.7);

    //--- define some media
    TGeoMedium *vacuum  = new TGeoMedium("Vacuum", 1, mat_vacuum);
    TGeoMedium *Al      = new TGeoMedium("Al", 2, mat_Al);

    //create the world box
    TGeoVolume *top = fGeom->MakeBox("top", vacuum, 200., 200., 200.);     
    top->SetVisibility(false); 
    fGeom->SetTopVolume(top); 

    //now, get ready to make a box for the sieve
    const double sieve_depth    = 12.7  /2.; 
    const double sieve_width    = 107.95/2.; 
    const double sieve_height   = 82.55 /2.;

    TGeoVolume *sieve_vol = fGeom->MakeBox("sieve_box", vacuum, sieve_width, sieve_height, sieve_depth); 
    sieve_vol->SetVisibility(kFALSE); 

    TGeoVolume *sieve = fGeom->MakeBox("sieve", Al, sieve_width, sieve_height, sieve_depth);     

    //two possible hole radii
    const double holeR_small = 0.6985; 
    const double holeR_big   = 1.3462; 

    TGeoVolume *hole_small = fGeom->MakeTube("small_hole", vacuum, holeR_small, holeR_small, 0.);
    //hole_small->SetLineColor(kBlack); 
    //hole_small->SetVisibility(kTRUE); 

    TGeoVolume *hole_big   = fGeom->MakeTube("big_hole",   vacuum, holeR_big, holeR_big, 0.); 

    //now, we're ready to draw the sieve-holes. 
    vector<SieveHole> sieve_holes = ApexOptics::ConstructSieveHoles(is_RHRS); 

    int i_hole_small=1;
    int i_hole_big=1; 

    for (SieveHole hole : sieve_holes) {   
        
        double x = hole.x * 1e3; 
        double y = hole.y * 1e3; 

        if (hole.is_big) {
            sieve_vol->AddNodeOverlap(hole_big,   i_hole_big++,   new TGeoTranslation(x,y, -sieve_depth)); 
        } else {
            sieve_vol->AddNodeOverlap(hole_small, i_hole_small++, new TGeoTranslation(x,y, -sieve_depth)); 
        }
    }
    sieve_vol->AddNodeOverlap(sieve, 1); 

    //create a set of axes 
    const double axes_size = 50.; 
    const unsigned int axis_color = kBlack; 
    TGeoVolume *xyz_axes = fGeom->MakeBox("axes", vacuum, axes_size, axes_size, axes_size); 
    xyz_axes->SetVisibility(kFALSE); 
    {
        TGeoVolume *single_axis = fGeom->MakeBox("single_axis", vacuum, axes_size, axes_size/10., axes_size/10.); 
        single_axis->SetVisibility(kFALSE); 
        
        TGeoRotation *single_axis_rot = new TGeoRotation; 
        single_axis_rot->SetAngles(90., 90., 0.); 

        TGeoVolume *axis_stick = fGeom->MakeTube("axis_stick", Al, axes_size/100., axes_size/100., axes_size*0.8); 
        axis_stick->SetLineColor(axis_color); 

        TGeoVolume *axis_cone  = fGeom->MakeCone("axis_cone", Al, axes_size*0.2, 0., axes_size*0.1, 0., 0.);
        axis_cone ->SetLineColor(axis_color);  
        
        single_axis->AddNodeOverlap(axis_stick, 1, new TGeoCombiTrans(axes_size*0.8, 0., 0., single_axis_rot)); 
        single_axis->AddNodeOverlap(axis_cone,  1, new TGeoCombiTrans(axes_size*1.8, 0., 0., single_axis_rot)); 

        auto rot_xaxis = new TGeoRotation; rot_xaxis->SetAngles(0., 0., 0.); 
        auto rot_yaxis = new TGeoRotation; rot_yaxis->SetAngles(90., 0., 0.); 
        auto rot_zaxis = new TGeoRotation; rot_zaxis->SetAngles(0., 90., 90.);
        
        //create letters
        const double letter_size = axes_size/10.; 
        //long letter stem
        TGeoVolume *letter_barL = fGeom->MakeBox("xbarL", Al, letter_size, axes_size/100., axes_size/100.); 
        letter_barL->SetLineColor(axis_color);
        //short letter stem
        TGeoVolume *letter_barS = fGeom->MakeBox("xbarS", Al, letter_size/2., axes_size/100., axes_size/100.); 
        letter_barS->SetLineColor(axis_color);

        //rotations to be used for letters
        auto rot_45p = new TGeoRotation; rot_45p->SetAngles(45., 0., 0.); 
        auto rot_0   = new TGeoRotation; rot_0  ->SetAngles(0., 0., 0.); 
        auto rot_45m = new TGeoRotation; rot_45m->SetAngles(-45., 0., 0.);
        auto rot_90  = new TGeoRotation; rot_90 ->SetAngles(90., 0., 0.); 

        //create 'X' 
        TGeoVolume *name_x = fGeom->MakeBox("X", vacuum, letter_size, letter_size, letter_size);
        name_x->SetVisibility(kFALSE); 
         
        name_x->AddNodeOverlap(letter_barL, 1, rot_45m); 
        name_x->AddNodeOverlap(letter_barL, 2, rot_45p); 

        //create 'Y' 
        TGeoVolume *name_y = fGeom->MakeBox("Y", vacuum, letter_size, letter_size, letter_size); 
        name_y->SetVisibility(kFALSE); 

        name_y->AddNodeOverlap(letter_barS, 1, new TGeoCombiTrans(  letter_size/(2.*sqrt(2.)), -letter_size/(2.*sqrt(2.)), 0., rot_45m)); 
        name_y->AddNodeOverlap(letter_barS, 2, new TGeoCombiTrans(  letter_size/(2.*sqrt(2.)), +letter_size/(2.*sqrt(2.)), 0., rot_45p)); 
        name_y->AddNodeOverlap(letter_barS, 2, new TGeoCombiTrans( -letter_size/2., 0., 0., rot_0)); 

        //create 'Z'
        TGeoVolume *name_z = fGeom->MakeBox("Z", vacuum, letter_size, letter_size, letter_size); 
        name_z->SetVisibility(kFALSE); 

        name_z->AddNodeOverlap(letter_barS, 1, new TGeoCombiTrans( +letter_size/2., 0., 0., rot_90)); 
        name_z->AddNodeOverlap(letter_barS, 2, new TGeoCombiTrans( -letter_size/2., 0., 0., rot_90)); 
        auto letter_barLsqrt2 = fGeom->MakeBox("xbarLsqrt2", Al, letter_size*sqrt(2.), axes_size/100., axes_size/100.);
        name_z->AddNodeOverlap(letter_barL, 1, new TGeoCombiTrans( 0., 0., 0., rot_45p)); 

        xyz_axes->AddNodeOverlap(single_axis, 1, rot_xaxis); xyz_axes->AddNodeOverlap(name_x, 1, new TGeoCombiTrans(axes_size*2.2, 0., 0., rot_xaxis)); 
        xyz_axes->AddNodeOverlap(single_axis, 2, rot_yaxis); xyz_axes->AddNodeOverlap(name_y, 1, new TGeoCombiTrans(0., axes_size*2.2, 0., rot_yaxis)); 
        xyz_axes->AddNodeOverlap(single_axis, 3, rot_zaxis); xyz_axes->AddNodeOverlap(name_z, 1, new TGeoCombiTrans(0., 0., axes_size*2.2, rot_zaxis)); 
    }
    //5trrrrrr610200000000000
    // -muon 

    //now, we're ready to add the sieve-holes
    
    //this fcn returns the value to us in radians, so we need to convert to degrees. 
    const double sieve_angle = ApexOptics::Get_sieve_angle(is_RHRS) * (180./TMath::Pi()); 
    
    TGeoRotation *rot_sieve = new TGeoRotation; 
    rot_sieve->SetAngles(90. + sieve_angle, 90., -90.); 
    
    //convert form meters to mm
    auto sieve_pos = ApexOptics::Get_sieve_pos(is_RHRS); 

    //set the z-coordinate to 0, before we convert to HCS
    sieve_pos[2] = 0.; 
    sieve_pos = ApexOptics::SCS_to_HCS(is_RHRS, sieve_pos) * 1e3; 

    auto trans_sieve = new TGeoCombiTrans( sieve_pos.z(), sieve_pos.x(), sieve_pos.y(), rot_sieve); 
    
    top->AddNodeOverlap(sieve_vol, 1, trans_sieve); 

    //draw HCS coordinate axes
    auto rot_HCS = new TGeoRotation; 
    rot_HCS->SetAngles(90., 90., 0.); 

    top->AddNodeOverlap(xyz_axes, 1, rot_HCS); 



    fGeom->CloseGeometry();  

    
    
    //now, we're ready to draw the geometry in our interactive window. 
    fCanvFrame = new TGVerticalFrame(this, 1600, 800); 
    
    fECanvas = new TRootEmbeddedCanvas("canvas_geom", fCanvFrame, 1600, 800); 

    fCanvFrame->AddFrame(fECanvas, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsCenterY | kLHintsExpandY, 0, 0, 5, 5)); 

    auto canv = fECanvas->GetCanvas();
    
    canv->cd(); 
    top->Draw(); 

    AddFrame(fCanvFrame, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsCenterY | kLHintsExpandY, 0, 0, 0, 0)); 

    SetWindowName("Sieve geometry");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
}
//_______________________________________________________________________________________________________________________________

int test_3d_geom(const bool is_RHRS=false)
{
    new GeometryFrame(gClient->GetRoot(), 1600, 800, is_RHRS); 
    return 0; 
}