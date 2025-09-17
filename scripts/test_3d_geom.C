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
#include <stdexcept> 
#include <iostream> 
#include <algorithm> 
#include <cstdio>
#include <TVirtualGeoTrack.h> 
#include <TGeoTrack.h> 
#include <TParticle.h> 
#include <ROOT/RVec.hxx>

//TGMainFrame windows
#include <TGWindow.h> 
#include <TGFrame.h> 
#include <TRootEmbeddedCanvas.h>
#include <TGButton.h>
#include <TGClient.h>  

using namespace std; 
using ApexOptics::OpticsTarget_t; 
using ApexOptics::Trajectory_t; 
using namespace ROOT::VecOps; 


//this is an enum to track which objects are drawn 
struct Object_t { 
    enum NameBit : char16_t {
        kNone       = 0, 
        kSieve_LHRS = 1 << 0,           //left sieve 
        kSieve_RHRS = 1 << 1,           //right sieve
        kBeam       = 1 << 2,           //beam 
        kTarget_production = 1 << 3,    //production foils
        kTarget_VWires     = 1 << 4,
        kTrack_L_real      = 1 << 5,    // LHRS - actual track
        kTrack_L_fg        = 1 << 6,    // LHRS - first-guess reconstruction
        kTrack_L_fan       = 1 << 7,    // LHRS - spread of possible track trajectories
        kHCS_axes          = 1 << 8
    };
    
    NameBit name; 
    TGeoVolume* volume=nullptr; 

    bool operator==(const Object_t& rhs) const { return rhs.name == name; }
};
//the definition of our geometry class
//_______________________________________________________________________________________________________________________________
class GeometryFrame : public TGMainFrame {
private: 
    TGHorizontalFrame *fCanvFrame; 
    TGHorizontalFrame *fButtonFrame; 

    TRootEmbeddedCanvas *fECanvas; 

    TGeoManager *fGeom; 
    TGeoVolume *fTopVolume; 

    TGeoMedium *fMedVacuum; 
    TGeoMedium *fMedAluminum; 

    TGeoRotation *fRot_HCS; 

    Object_t::NameBit fDrawnObjects = Object_t::kNone; 

    std::vector<Object_t> fObjects{}; 
    std::vector<std::pair<Object_t::NameBit,TGCheckButton*>> fButtons{}; 

    //a map of buttons, whith their associated name-bit. 
    Object_t* FindObject(Object_t::NameBit name); 

    //create geometry
    void CreateGeometry(); 

    //create a set of axes
    TGeoVolume* CreateAxesGeometry(const double axes_size=50., unsigned int color=kBlack); 

    //create geometry for the left (or right) sieve
    void CreateSieveGeometry(const bool arm_is_RHRS); 
    
    //create geometry for the production target
    void CreateProductionTargetGeometry(); 

    //create geometry for the VWires 
    void CreateVWireGeometry(unsigned int color=kBlack); 

    //create geometry for the beam
    void CreateBeamGeometry(const double x_hcs, const double y_hcs, unsigned int color=kRed); 

    //create track geometry
    struct Track_t { 
        double x_sv, y_sv; 
        TVector3 vtx_hcs; 
        bool is_RHRS=false; 
        TGeoTrack* geo_track=nullptr; 
    }; 
    void CreateTrackGeometry(const bool is_RHRS, GeometryFrame::Track_t& track, unsigned int color=kRed);

    std::vector<Track_t> fTracks_fan{}; 
    Track_t fTrack_real; 
    Track_t fTrack_first_guess; 

    const bool is_RHRS; 

public: 
    GeometryFrame(const TGWindow* p, UInt_t w, UInt_t h, const bool is_RHRS=false); 
    ~GeometryFrame(); 

    void ButtonClicked(); 

    ClassDef(GeometryFrame,1); 
};
//_______________________________________________________________________________________________________________________________

//_______________________________________________________________________________________________________________________________
Object_t* GeometryFrame::FindObject(Object_t::NameBit name)
{
    //find an object with the given namebit 
    auto it = std::find_if( fObjects.begin(), fObjects.end(), [name](const Object_t& obj){ return obj.name == name; });

    if (it == fObjects.end()) {
        throw invalid_argument("in <GeometryFrame::FindObject>: invalid object name bit given."); 
        return nullptr; 
    }

    return &(*it); 
}
//_______________________________________________________________________________________________________________________________
void GeometryFrame::CreateBeamGeometry(const double x_hcs, const double y_hcs, unsigned int color=kRed)
{
    //draw the beam
    const double beam_length = 500.; 
    const double beam_rad    = 0.10; 

    const double cone_rad    = 2.; 

    TGeoVolume *vol = fGeom->MakeBox("beam_box", fMedVacuum, 20., 20., beam_length/2.);
    vol->SetVisibility(kFALSE);

    TGeoVolume *beam_stick = fGeom->MakeTube("beam_stick", fMedVacuum, beam_rad, beam_rad, beam_length); 
    TGeoVolume *beam_cone  = fGeom->MakeCone("beam_cone",  fMedVacuum, cone_rad, 0., cone_rad, 0., 0.); 

    beam_stick->SetLineColor(color); 
    beam_cone ->SetLineColor(color); 

    vol->AddNodeOverlap(beam_stick, 1); 
    vol->AddNodeOverlap(beam_cone,  1, new TGeoTranslation(0., 0., beam_length)); 

    fTopVolume->AddNodeOverlap(vol, 1, new TGeoCombiTrans(250., 0., 0., fRot_HCS)); 

    fObjects.push_back({.name=Object_t::kBeam, .volume=vol}); 

    return; 
}
//_______________________________________________________________________________________________________________________________
void GeometryFrame::CreateVWireGeometry(const unsigned int color)
{
    //create geometry for the v-wire targets. 
    //recall that all units are in mm
    OpticsTarget_t V1 = ApexOptics::GetTarget("V1"); 
    OpticsTarget_t V2 = ApexOptics::GetTarget("V2"); 
    OpticsTarget_t V3 = ApexOptics::GetTarget("V3"); 
    
    const double wire_radius = 0.050; 
    const double wire_height = 50.; 

    TGeoVolume *vol = fGeom->MakeBox("vwire_box", fMedVacuum, 5., 50., 0.5*1e3*(V3.z_hcs - V1.z_hcs)); 
    vol->SetVisibility(kFALSE); 

    TGeoVolume *wire = fGeom->MakeTube("vwire", fMedAluminum, wire_radius, wire_radius, wire_height/2.); 
    wire->SetLineColor(color); 

    auto rot_zaxis = new TGeoRotation; rot_zaxis->SetAngles(0., 90., 90.);

    vol->AddNodeOverlap(wire, 1, new TGeoCombiTrans(1e3*V1.x_hcs, 0., 1e3*V1.z_hcs, rot_zaxis)); 
    vol->AddNodeOverlap(wire, 2, new TGeoCombiTrans(1e3*V2.x_hcs, 0., 1e3*V2.z_hcs, rot_zaxis)); 
    vol->AddNodeOverlap(wire, 3, new TGeoCombiTrans(1e3*V3.x_hcs, 0., 1e3*V3.z_hcs, rot_zaxis)); 

    fTopVolume->AddNodeOverlap(vol, 1, fRot_HCS); 

    fObjects.push_back({.name=Object_t::kTarget_VWires, .volume=vol}); 

    return; 
}
//_______________________________________________________________________________________________________________________________
void GeometryFrame::CreateProductionTargetGeometry()
{
    //create geometry for the 'produciton' target.
    //recall that all units are in mm
    const int n_targets = 10.; 
    const double target_width = 2.5; 
    const double target_depth = 0.010; 
    const double target_spacing = 55.; 

    const double z_foil_first = -238.19; 
    
    //the targets are much taller than this, but there's not much point in drawing that, as the produciton raster amplitude is ~2mm
    const double target_height = 50.; 

    TGeoVolume *prod_vol = fGeom->MakeBox("prod_box", fMedVacuum, target_width/2., target_height/2., 0.5*(target_spacing * (double)n_targets)); 
    prod_vol->SetVisibility(kFALSE); 

    int i_foil=1; 

    const unsigned int color_prod = kBlack; 

    double z_targ = z_foil_first; 
    
    auto target = fGeom->MakeBox("prod_foil", fMedAluminum, target_width/2., target_height/2., target_depth/2.); 
    target->SetLineColor(color_prod); 

    for (unsigned int i_targ=1; i_targ<=n_targets; i_targ++) {

        prod_vol->AddNodeOverlap(target, i_targ, new TGeoTranslation(0., 0., z_targ)); 

        z_targ += target_spacing; 
    }

    fTopVolume->AddNodeOverlap(prod_vol, 1, fRot_HCS); 

    fObjects.push_back({.name=Object_t::kTarget_production, .volume=prod_vol}); 

    return; 
} 
//_______________________________________________________________________________________________________________________________
void GeometryFrame::CreateSieveGeometry(const bool arm_is_RHRS)
{   
    //now, get ready to make a box for the sieve
    const double sieve_depth    = 12.7  /2.; 
    const double sieve_width    = 107.95/2.; 
    const double sieve_height   = 82.55 /2.;

    TGeoVolume *sieve_vol = fGeom->MakeBox(Form("sieve_box_%s",arm_is_RHRS?"R":"L"), fMedVacuum, sieve_width, sieve_height, sieve_depth); 
    sieve_vol->SetVisibility(kFALSE); 

    TGeoVolume *sieve = fGeom->MakeBox("sieve", fMedAluminum, sieve_width, sieve_height, sieve_depth);     

    //two possible hole radii
    const double holeR_small = 0.6985; 
    const double holeR_big   = 1.3462; 

    TGeoVolume *hole_small = fGeom->MakeTube("small_hole", fMedVacuum, holeR_small, holeR_small, 0.);
    //hole_small->SetLineColor(kBlack); 
    //hole_small->SetVisibility(kTRUE); 

    TGeoVolume *hole_big   = fGeom->MakeTube("big_hole",   fMedVacuum, holeR_big, holeR_big, 0.); 

    //now, we're ready to draw the sieve-holes. 
    vector<SieveHole> sieve_holes = ApexOptics::ConstructSieveHoles(arm_is_RHRS); 

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

    Object_t::NameBit name = arm_is_RHRS ? Object_t::kSieve_RHRS : Object_t::kSieve_LHRS; 

    //this fcn returns the value to us in radians, so we need to convert to degrees. 
    const double sieve_angle = ApexOptics::Get_sieve_angle(arm_is_RHRS) * (180./TMath::Pi()); 
    
    TGeoRotation *rot_sieve = new TGeoRotation; 
    rot_sieve->SetAngles(90. + sieve_angle, 90., -90.); 
    
    //convert form meters to mm
    auto sieve_pos = ApexOptics::Get_sieve_pos(arm_is_RHRS); 

    //set the z-coordinate to 0, before we convert to HCS
    sieve_pos[2] = 0.; 
    sieve_pos = ApexOptics::SCS_to_HCS(arm_is_RHRS, sieve_pos) * 1e3; 

    auto trans_sieve = new TGeoCombiTrans( sieve_pos.z(), sieve_pos.x(), sieve_pos.y(), rot_sieve); 
    
    fTopVolume->AddNodeOverlap(sieve_vol, 1, trans_sieve); 

    fObjects.push_back({.name=name, .volume=sieve_vol}); 

    return; 
}
//_______________________________________________________________________________________________________________________________
void GeometryFrame::ButtonClicked()
{
    fDrawnObjects = Object_t::kNone; 
    for (auto& button : fButtons) {
        if (button.second && button.second->IsOn()) fDrawnObjects = Object_t::NameBit( fDrawnObjects | button.first ); 
    }

    for (auto & obj : fObjects) {
        if (!obj.volume) continue; 
        //if the bit for this object is set to 'true', then draw it. 
        if (fDrawnObjects & obj.name) { 
            obj.volume->SetVisDaughters(kTRUE);
        } else { 
            obj.volume->SetVisDaughters(kFALSE); 
        }
    }

    if (fECanvas->GetCanvas()) fECanvas->GetCanvas()->cd(); 
    if (fTopVolume) fTopVolume->Draw();   

    const unsigned int color_real   = kRed; 
    const unsigned int color_fg     = kBlue; 
    const unsigned int color_fan    = kBlue; 

    //handle track drawing
    // -- real
    if (fDrawnObjects & Object_t::kTrack_L_real) {
        fTrack_real.geo_track->SetLineColorAlpha(color_real, 1.); 
        fTrack_real.geo_track->Draw(); 
    }

    // -- fan
    if (fDrawnObjects & Object_t::kTrack_L_fan) {
        for (auto& track : fTracks_fan) {
            track.geo_track->SetLineColorAlpha(color_real, 1.); 
            track.geo_track->Draw(); 
        }
    }

    // -- first-guess forward
    if (fDrawnObjects & Object_t::kTrack_L_fg) {
        fTrack_first_guess.geo_track->SetLineColorAlpha(color_fan, 0.15); 
        fTrack_first_guess.geo_track->Draw();
    } 
    
    return; 
}
//_______________________________________________________________________________________________________________________________

TGeoVolume* GeometryFrame::CreateAxesGeometry(const double axes_size=50., const unsigned int color)
{
    //create a set of axes, with the given orientation, position, size and color. 
    //create a set of axes 
    TGeoVolume *xyz_axes = fGeom->MakeBox("axes", fMedVacuum, axes_size, axes_size, axes_size); 
    xyz_axes->SetVisibility(kFALSE); 
    
    TGeoVolume *single_axis = fGeom->MakeBox("single_axis", fMedVacuum, axes_size, axes_size/10., axes_size/10.); 
    single_axis->SetVisibility(kFALSE); 
    
    TGeoRotation *single_axis_rot = new TGeoRotation; 
    single_axis_rot->SetAngles(90., 90., 0.); 

    TGeoVolume *axis_stick = fGeom->MakeTube("axis_stick", fMedAluminum, axes_size/100., axes_size/100., axes_size*0.8); 
    axis_stick->SetLineColor(color); 

    TGeoVolume *axis_cone  = fGeom->MakeCone("axis_cone", fMedAluminum, axes_size*0.2, 0., axes_size*0.1, 0., 0.);
    axis_cone ->SetLineColor(color);  
    
    single_axis->AddNodeOverlap(axis_stick, 1, new TGeoCombiTrans(axes_size*0.8, 0., 0., single_axis_rot)); 
    single_axis->AddNodeOverlap(axis_cone,  1, new TGeoCombiTrans(axes_size*1.8, 0., 0., single_axis_rot)); 

    auto rot_xaxis = new TGeoRotation; rot_xaxis->SetAngles(0., 0., 0.); 
    auto rot_yaxis = new TGeoRotation; rot_yaxis->SetAngles(90., 0., 0.); 
    auto rot_zaxis = new TGeoRotation; rot_zaxis->SetAngles(0., 90., 90.);
    
    //create letters
    const double letter_size = axes_size/10.; 
    //long letter stem
    TGeoVolume *letter_barL = fGeom->MakeBox("xbarL", fMedAluminum, letter_size, axes_size/100., axes_size/100.); 
    letter_barL->SetLineColor(color);
    //short letter stem
    TGeoVolume *letter_barS = fGeom->MakeBox("xbarS", fMedAluminum, letter_size/2., axes_size/100., axes_size/100.); 
    letter_barS->SetLineColor(color);

    //rotations to be used for letters
    auto rot_45p = new TGeoRotation; rot_45p->SetAngles(45., 0., 0.); 
    auto rot_0   = new TGeoRotation; rot_0  ->SetAngles(0., 0., 0.); 
    auto rot_45m = new TGeoRotation; rot_45m->SetAngles(-45., 0., 0.);
    auto rot_90  = new TGeoRotation; rot_90 ->SetAngles(90., 0., 0.); 

    //create 'X' 
    TGeoVolume *name_x = fGeom->MakeBox("X", fMedVacuum, letter_size, letter_size, letter_size);
    name_x->SetVisibility(kFALSE); 
        
    name_x->AddNodeOverlap(letter_barL, 1, rot_45m); 
    name_x->AddNodeOverlap(letter_barL, 2, rot_45p); 

    //create 'Y' 
    TGeoVolume *name_y = fGeom->MakeBox("Y", fMedVacuum, letter_size, letter_size, letter_size); 
    name_y->SetVisibility(kFALSE); 

    name_y->AddNodeOverlap(letter_barS, 1, new TGeoCombiTrans(  letter_size/(2.*sqrt(2.)), -letter_size/(2.*sqrt(2.)), 0., rot_45m)); 
    name_y->AddNodeOverlap(letter_barS, 2, new TGeoCombiTrans(  letter_size/(2.*sqrt(2.)), +letter_size/(2.*sqrt(2.)), 0., rot_45p)); 
    name_y->AddNodeOverlap(letter_barS, 2, new TGeoCombiTrans( -letter_size/2., 0., 0., rot_0)); 

    //create 'Z'
    TGeoVolume *name_z = fGeom->MakeBox("Z", fMedVacuum, letter_size, letter_size, letter_size); 
    name_z->SetVisibility(kFALSE); 

    name_z->AddNodeOverlap(letter_barS, 1, new TGeoCombiTrans( +letter_size/2., 0., 0., rot_90)); 
    name_z->AddNodeOverlap(letter_barS, 2, new TGeoCombiTrans( -letter_size/2., 0., 0., rot_90)); 
    name_z->AddNodeOverlap(letter_barL, 1, new TGeoCombiTrans( 0., 0., 0., rot_45p)); 

    xyz_axes->AddNodeOverlap(single_axis, 1, rot_xaxis); xyz_axes->AddNodeOverlap(name_x, 1, new TGeoCombiTrans(axes_size*2.2, 0., 0., rot_xaxis)); 
    xyz_axes->AddNodeOverlap(single_axis, 2, rot_yaxis); xyz_axes->AddNodeOverlap(name_y, 1, new TGeoCombiTrans(0., axes_size*2.2, 0., rot_yaxis)); 
    xyz_axes->AddNodeOverlap(single_axis, 3, rot_zaxis); xyz_axes->AddNodeOverlap(name_z, 1, new TGeoCombiTrans(0., 0., axes_size*2.2, rot_zaxis)); 

    return xyz_axes; 
}
//_______________________________________________________________________________________________________________________________
void GeometryFrame::CreateTrackGeometry(const bool is_RHRS, GeometryFrame::Track_t& real_track, unsigned int color)
{
    const char* const here = "draw_trajectory_fan"; 

    const vector<string> branches_sv{
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    }; 

    const vector<string> branches_fp{
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    const char* path_poly_fp_sv = "data/csv/poly_prod_fp_sv_L_6ord.dat"; 
    const char* path_poly_sv_fp = "data/csv/poly_prod_sv_fp_L_6ord.dat"; 

    //try to parse the NPolyArray-s from the given db-files
    NPolyArray parr_fp_sv, parr_sv_fp; 
    try {
        parr_fp_sv = ApexOptics::Parse_NPolyArray_from_file(path_poly_fp_sv, branches_sv, 4); 
        parr_sv_fp = ApexOptics::Parse_NPolyArray_from_file(path_poly_sv_fp, branches_fp, 5); 
    
    } catch (const std::exception& e) {

        Error(here, "Exception caught trying to parse NPolyArray db-files.\n what(): %s", e.what()); 
        return; 
    }

    //this function takes a starting vertex (in HCS), and a target position on the sieve-face, and given dp/p momentum parameter, 
    //and generates a 'Trajectory_t' in SCS 
    auto Generate_test_trajectory = [&](TVector3 react_vertex, double x_sv, double y_sv, double dpp)
    {
        //first, convert the TVector3 to SCS. 
        react_vertex.RotateY( -ApexOptics::Get_sieve_angle(is_RHRS) ); 
        react_vertex.RotateZ( TMath::Pi()/2. ); 

        react_vertex += -ApexOptics::Get_sieve_pos(is_RHRS);

        //now compute the 'angles' in scs, using the intercept with the sieve-face
        double dxdz = ( x_sv - react_vertex.x() )/( 0. - react_vertex.z() ); 
        double dydz = ( y_sv - react_vertex.y() )/( 0. - react_vertex.z() ); 

        //return a Trajectory_t struct with all this information 
        return Trajectory_t{ 
            .x    = x_sv, 
            .y    = y_sv, 
            .dxdz = dxdz, 
            .dydz = dydz, 
            .dpp  = dpp
        }; 
    }; 

    

    const double center_dpp = +0.000; 

    vector<Trajectory_t> test_trajectories{};

    TVector3 react_vtx_center = real_track.vtx_hcs * 1e-3; 

    const double z_target = react_vtx_center.z(); 


    Trajectory_t real_traj = Generate_test_trajectory(react_vtx_center, real_track.x_sv, real_track.y_sv, center_dpp); 

    //the max +/- level of dp/p to search. 
    const double dpp_search_range = 0.50e-3; 

    //double this number +1 is the nubmer of 'trajectory points' to draw
    const size_t half_trajectory_points = 20; 

    //first, try to find the 'focal plane coordinate' cooresponding to this trajectory. 

    auto Xfp_rvec = parr_sv_fp.Eval(ApexOptics::Trajectory_t_to_RVec(real_traj));  

    auto Xsv_fwd_rvec = parr_fp_sv.Eval(Xfp_rvec); 

    const RVec<double> Xsv_rvec_actual = ApexOptics::Trajectory_t_to_RVec(real_traj); 

    Trajectory_t traj_fwd_model = ApexOptics::RVec_to_Trajectory_t(Xsv_fwd_rvec);

    //cout << "new traj.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
    //printf("traj (real vs fwd-model:)  {x,y,dxdz,dydz,dpp}:\n"); 
    //printf(" -- {%+.4f,%+.4f,%+.4f,%+.4f,%+.4f}\n", traj.x, traj.y, traj.dxdz, traj.dydz, traj.dpp); 
    //printf(" -- {%+.4f,%+.4f,%+.4f,%+.4f,%+.4f}\n", traj_fwd_model.x, traj_fwd_model.y, traj_fwd_model.dxdz, traj_fwd_model.dydz, traj_fwd_model.dpp); 
    
    const double error_threshold = 1e-9; 

    auto rvec_mag = [](const RVec<double>& V) { double val=0.; for(const auto& x : V) { val += x*x; } return sqrt(val); }; 

    //now, iterate to a Xsv coordinate which 'matches' the Xfp coordinate. 
    auto Xsv_rvec = Xsv_fwd_rvec;
    
    printf(" after %i iterations, final error: % .3e\n", 
        parr_sv_fp.Iterate_to_root(Xsv_rvec, Xfp_rvec, 15, error_threshold), 
        rvec_mag(Xfp_rvec - parr_sv_fp.Eval(Xsv_rvec))
    ); 

    //now, let's do our little dance of finding 'adjacent' trajectories to draw. 
    auto Find_next_traj = [&](RVec<double> Xsv, const double d_dpp) 
    {
        auto J = parr_sv_fp.Jacobian(Xsv_rvec); 

        //the last column of this jacobian corresponds to dp/p. since this is the variable we will be varying, we will 
        RVec<double> J4{ 
            -J.at(0,4),
            -J.at(1,4), 
            -J.at(2,4), 
            -J.at(3,4)
        }; 

        RMatrix Ji(4,4, {
            J.at(0,0), J.at(0,1), J.at(0,2), J.at(0,3),
            J.at(1,0), J.at(1,1), J.at(1,2), J.at(1,3),
            J.at(2,0), J.at(2,1), J.at(2,2), J.at(2,3),
            J.at(3,0), J.at(3,1), J.at(3,2), J.at(3,3)
        }); 

        auto dX = Ji.Solve( J4 ); 

        if (dX.empty()) return RVec<double>{}; 
        for (double x : dX) { if (x != x) return RVec<double>{}; }
        
        dX *= d_dpp; 

        dX.push_back( d_dpp ); 

        Xsv += dX; 

        //now, iterate to 'correct' whatever small error may be in our linear extrapolation 'Xsv += dX' 
        parr_sv_fp.Iterate_to_root(Xsv, Xfp_rvec, 15, 1e-10); 

        return Xsv; 
    }; 
    //_________________________________________________________________________________________________
    
    
    const double d_dpp = dpp_search_range/((double)half_trajectory_points-1); 

    //if the next vector is further than this, then quit finding new points (something is wrong!)
    const double max_dist_to_next = d_dpp * 20.; 

    vector<Trajectory_t> trajectories_forward; 

    auto Xsv_next = Xsv_rvec; 
    
    for (size_t i=0; i<half_trajectory_points; i++) {
        
        Xsv_next = Find_next_traj(Xsv_next, d_dpp); 

        //this happens if the 'Find_next_traj' fails for some reason. if so, stop finding new points. 
        if (Xsv_next.size() != 5) break; 

        trajectories_forward.push_back( ApexOptics::RVec_to_Trajectory_t(Xsv_next) ); 
    }

    //now, find the 'backward' trajectories
    vector<Trajectory_t> trajectories_backward;
    
    Xsv_next = Xsv_rvec; 

    for (size_t i=0; i<half_trajectory_points; i++) {

        auto Xsv_new = Find_next_traj(Xsv_next, -1.*d_dpp); 

        //this happens if the 'Find_next_traj' fails for some reason. if so, stop finding new points. 
        if (Xsv_next.size() != 5) break; 
        
        //check if the new vector has changed too much. if so, something has gone wrong in the iteration process
        if (rvec_mag(Xsv_next - Xsv_new) > max_dist_to_next) break; 

        //check if the momentum is 'out of range' of the maximum search range
        if (fabs(Xsv_next[4] - real_traj.dpp) > dpp_search_range) break; 

        Xsv_next = Xsv_new; 

        trajectories_backward.push_back( ApexOptics::RVec_to_Trajectory_t(Xsv_next) );
    }

    //now we're done. make the vectors. 
    cout << "number of trajectories (fwd/back): " 
            << trajectories_forward.size() << " / " 
            << trajectories_backward.size() << endl; 

    //now, combine them all in one vector. 
    vector<Trajectory_t> trajectories;
    trajectories.reserve( trajectories_forward.size() + trajectories_backward.size() + 1 );
    
    for (int i=trajectories_backward.size()-1; i>=0; i--) 
        trajectories.push_back( trajectories_backward.at(i) ); 
    
    trajectories.push_back( ApexOptics::RVec_to_Trajectory_t(Xsv_rvec) ); 

    for (int i=0; i<trajectories_forward.size(); i++) 
        trajectories.push_back( trajectories_forward.at(i) ); 

    
    //now, convert this data to tracks. 
    int i_track=100; 

    auto Trajectory_t_to_Track_t = [&](Trajectory_t traj)
    {
        Track_t track; 

        track.x_sv = traj.x * 1e3; 
        track.y_sv = traj.y * 1e3; 

        TVector3 pt_sieve = ApexOptics::SCS_to_HCS(is_RHRS, TVector3(traj.x, traj.y, 0.)); 

        auto traj_hcs = ApexOptics::SCS_to_HCS(is_RHRS, traj); 

        double x_target = traj_hcs.x + (traj_hcs.dxdz * z_target); 
        double y_target = traj_hcs.y + (traj_hcs.dydz * z_target);     

        track.vtx_hcs = TVector3(x_target, y_target, z_target) * 1e3; 
        track.is_RHRS = is_RHRS; 

        track.geo_track = (TGeoTrack*)fGeom->MakeTrack(i_track++, 0., nullptr); 

        track.geo_track->AddPoint(track.vtx_hcs.z(), track.vtx_hcs.x(), track.vtx_hcs.y(), 0.); 
        track.geo_track->AddPoint(1e3*pt_sieve.z(),  1e3*pt_sieve.x(),  1e3*pt_sieve.y(),  1.); 

        return track; 
    };

    real_track         = Trajectory_t_to_Track_t(real_traj); 
    fTrack_first_guess = Trajectory_t_to_Track_t(traj_fwd_model); 

    fTracks_fan.clear(); 
    for (auto traj : trajectories) fTracks_fan.push_back( Trajectory_t_to_Track_t(traj) ); 

    fGeom->SetTminTmax(0., 1.); 
    
    //fGeom->DrawTracks(); 
    return; 
}
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________


GeometryFrame::~GeometryFrame() 
{
    fGeom->~TGeoManager(); 
    Cleanup(); 
}
//_______________________________________________________________________________________________________________________________
void GeometryFrame::CreateGeometry() 
{
    //noop  
    fGeom = new TGeoManager("simple1", "Simple geometry");

    //--- define some materials
    TGeoMaterial *mat_vacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *mat_Al     = new TGeoMaterial("Al", 26.98, 13, 2.7);

    //--- define some media
    fMedVacuum   = new TGeoMedium("Vacuum", 1, mat_vacuum);
    fMedAluminum = new TGeoMedium("Al", 2, mat_Al);

    //create the world box
    fTopVolume = fGeom->MakeBox("top", fMedVacuum, 200., 200., 200.);     
    fTopVolume->SetVisibility(false); 
    fGeom->SetTopVolume(fTopVolume); 


    //5trrrrrr610200000000000
    // -muon     
    
    //draw HCS coordinate axes
    auto hcs_axes = CreateAxesGeometry(50., kBlack); 
    fRot_HCS = new TGeoRotation; 
    fRot_HCS->SetAngles(90., 90., 0.); 
    fTopVolume->AddNodeOverlap(hcs_axes, 1, new TGeoCombiTrans(0., 0., -200., fRot_HCS)); 
    
    //now, we're ready to add the sieve-holes
    CreateSieveGeometry(true);  //RHRS 
    CreateSieveGeometry(false); //LHRS

    //create the production target geometry 
    CreateProductionTargetGeometry(); 

    //create vwire geometry
    CreateVWireGeometry(kBlue); 
   
    //create the beam geometry
    CreateBeamGeometry(0., 0.); 

    //create the track geometry
    fTrack_real = Track_t{.x_sv=-25e-3, .y_sv=-25e-3, .vtx_hcs=TVector3(0., 0., 36.281)}; 
    CreateTrackGeometry(is_RHRS, fTrack_real); 
    
    fGeom->CloseGeometry();  
}
//_______________________________________________________________________________________________________________________________

//_______________________________________________________________________________________________________________________________
GeometryFrame::GeometryFrame(const TGWindow* p, UInt_t w, UInt_t h, const bool _is_RHRS)
    : TGMainFrame(p, w, h), is_RHRS(_is_RHRS)
{
    CreateGeometry(); 

    //create the button frame
    fButtonFrame = new TGHorizontalFrame(this, 100, 800); 

    TGCheckButton *button=nullptr; 

    // -- LHRS sieve
    button = new TGCheckButton(this, "LHRS sieve");                             //add this button and name it
    button->Connect("Clicked()", "GeometryFrame", this, "ButtonClicked()");
    fButtonFrame->AddFrame(button, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
    fButtons.push_back({Object_t::kSieve_LHRS, button});                        //add it to the list of buttons
    button->SetState(kButtonDown);                                              //set its default state to 'on' 

    // -- RHRS sieve
    button = new TGCheckButton(this, "RHRS sieve");                             //add this button and name it
    button->Connect("Clicked()", "GeometryFrame", this, "ButtonClicked()");
    fButtonFrame->AddFrame(button, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
    fButtons.push_back({Object_t::kSieve_RHRS, button});                        //add it to the list of buttons
    button->SetState(kButtonDown);                                              //set its default state to 'on' 

    // -- production target
    button = new TGCheckButton(this, "Prod. Tungsten foils");                   //add this button and name it
    button->Connect("Clicked()", "GeometryFrame", this, "ButtonClicked()");
    fButtonFrame->AddFrame(button, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
    fButtons.push_back({Object_t::kTarget_production, button});                 //add it to the list of buttons
    button->SetState(kButtonDown);                                              //set its default state to 'on' 

    // -- vertical wire targets
    button = new TGCheckButton(this, "Vertical Optic Wires");                   //add this button and name it
    button->Connect("Clicked()", "GeometryFrame", this, "ButtonClicked()");
    fButtonFrame->AddFrame(button, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
    fButtons.push_back({Object_t::kTarget_VWires, button});                     //add it to the list of buttons
    button->SetState(kButtonUp);                                                //set its default state to 'on' 

    // -- beam
    button = new TGCheckButton(this, "Beam");                                   //add this button and name it
    button->Connect("Clicked()", "GeometryFrame", this, "ButtonClicked()");
    fButtonFrame->AddFrame(button, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
    fButtons.push_back({Object_t::kBeam, button});                              //add it to the list of buttons
    button->SetState(kButtonUp);                                                //set its default state to 'on' 

    // -- 'real' track
    button = new TGCheckButton(this, "Real Track");                             //add this button and name it
    button->Connect("Clicked()", "GeometryFrame", this, "ButtonClicked()");
    fButtonFrame->AddFrame(button, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
    fButtons.push_back({Object_t::kTrack_L_real, button});                      //add it to the list of buttons
    button->SetState(kButtonUp);                                                //set its default state to 'on' 

    // -- 'first guess' reconstructed track
    button = new TGCheckButton(this, "Forward Reco.");                          //add this button and name it
    button->Connect("Clicked()", "GeometryFrame", this, "ButtonClicked()");
    fButtonFrame->AddFrame(button, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
    fButtons.push_back({Object_t::kTrack_L_fg, button});                        //add it to the list of buttons
    button->SetState(kButtonUp);                                                //set its default state to 'on' 

    // -- 'fan' of reconstructed tracks
    button = new TGCheckButton(this, "Trajectory Spread");                      //add this button and name it
    button->Connect("Clicked()", "GeometryFrame", this, "ButtonClicked()");
    fButtonFrame->AddFrame(button, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
    fButtons.push_back({Object_t::kTrack_L_fan, button});                       //add it to the list of buttons
    button->SetState(kButtonUp);                                                //set its default state to 'on' 

    AddFrame(fButtonFrame, new TGLayoutHints(kLHintsLeft, 0,0,0,0)); 

    
    //now, we're ready to draw the geometry in our interactive window. 
    fCanvFrame = new TGHorizontalFrame(this, 1600, 800); 
    
    fECanvas = new TRootEmbeddedCanvas("canvas_geom", fCanvFrame, 1600, 800); 

    fCanvFrame->AddFrame(fECanvas, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY, 0, 0, 5, 5)); 

    auto canv = fECanvas->GetCanvas();
    
    canv->cd(); 
    //call this to draw only the default-drawn objects
    ButtonClicked(); 

    //this creates all relevant geometry

    AddFrame(fCanvFrame, new TGLayoutHints(kLHintsRight | kLHintsExpandX, 0, 0, 0, 0)); 

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