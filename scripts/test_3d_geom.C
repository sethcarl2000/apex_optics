#include <TGeoManager.h> 
#include <TGeoMaterial.h> 
#include <TGeoMedium.h> 
#include <TGeoMatrix.h> 
#include <TGeoArb8.h> 
#include <TSystem.h> 
#include <TEveManager.h> 
#include <TEvePointSet.h> 
#include <TEveArrow.h>
#include <TBrowser.h>  
#include <SieveHole.h> 
#include <ApexOptics.h>
#include <vector>  

using namespace std; 

int test_3d_geom(const bool is_RHRS=false)
{
    /*auto manager = new TGeoManager("world", "test"); 

    auto vacuum = new TGeoMaterial("Vacuum", 0,0,0); 

    auto medium_vac = new TGeoMedium("Vacuum",1,vacuum); 

    auto top = manager->MakeBox("Top", medium_vac, 10., 10., 10.); 

    manager->SetTopVolume(top); 

    manager->CloseGeometry(); 

    top->SetLineColor(kMagenta); 

    top->Draw(); 

    return 0;*/ 
    
    //gStyle->SetCanvasPreferGL(true)

    TGeoManager *geom = new TGeoManager("simple1", "Simple geometry");
    
    //--- define some materials
    TGeoMaterial *mat_vacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *mat_Al     = new TGeoMaterial("Al", 26.98, 13, 2.7);

    //--- define some media
    TGeoMedium *vacuum  = new TGeoMedium("Vacuum", 1, mat_vacuum);
    TGeoMedium *Al      = new TGeoMedium("Al", 2, mat_Al);

    //create the world box
    TGeoVolume *top = geom->MakeBox("top", vacuum, 200., 200., 200.);     
    top->SetVisibility(false); 
    geom->SetTopVolume(top); 

    

    //now, get ready to make a box for the sieve
    const double sieve_depth    = 12.7  /2.; 
    const double sieve_width    = 107.95/2.; 
    const double sieve_height   = 82.55 /2.;

    TGeoVolume *sieve_vol = geom->MakeBox("sieve_box", vacuum, sieve_width, sieve_height, sieve_depth); 
    sieve_vol->SetVisibility(kFALSE); 

    TGeoVolume *sieve = geom->MakeBox("sieve", Al, sieve_width, sieve_height, sieve_depth);     

    //two possible hole radii
    const double holeR_small = 0.6985; 
    const double holeR_big   = 1.3462; 

    TGeoVolume *hole_small = geom->MakeTube("small_hole", vacuum, holeR_small, holeR_small, 0.);
    //hole_small->SetLineColor(kBlack); 
    //hole_small->SetVisibility(kTRUE); 

    TGeoVolume *hole_big   = geom->MakeTube("big_hole",   vacuum, holeR_big, holeR_big, 0.); 

    //now, we're ready to draw the sieve-holes. 
    vector<SieveHole> sieve_holes = ApexOptics::ConstructSieveHoles(is_RHRS); 

    int i_hole_small=1;
    int i_hole_big=1; 

    for (SieveHole hole : sieve_holes) {   
        
        double x = hole.x * 1e3; 
        double y = hole.y * 1e3; 

        if (hole.is_big) {
            sieve_vol->AddNodeOverlap(hole_big,   i_hole_big++,   new TGeoTranslation(x,y, sieve_depth)); 
        } else {
            sieve_vol->AddNodeOverlap(hole_small, i_hole_small++, new TGeoTranslation(x,y, sieve_depth)); 
        }
    }
    sieve_vol->AddNodeOverlap(sieve, 1); 

    //create a set of axes 
    const double axes_size = 20.; 
    const unsigned int axis_color = kBlack; 
    TGeoVolume *xyz_axes = geom->MakeBox("axes", vacuum, axes_size, axes_size, axes_size); 
    xyz_axes->SetVisibility(kFALSE); 

    {
        TGeoVolume *single_axis = geom->MakeBox("single_axis", vacuum, axes_size, axes_size/10., axes_size/10.); 
        single_axis->SetVisibility(kFALSE); 
        
        TGeoRotation *single_axis_rot = new TGeoRotation; 
        single_axis_rot->SetAngles(90., 90., 0.); 

        TGeoVolume *axis_stick = geom->MakeTube("axis_stick", Al, axes_size/100., axes_size/100., axes_size*0.8); 
        axis_stick->SetLineColor(axis_color); 

        TGeoVolume *axis_cone  = geom->MakeCone("axis_cone", Al, axes_size*0.2, 0., axes_size*0.1, 0., 0.);
        axis_cone ->SetLineColor(axis_color);  
        
        single_axis->AddNodeOverlap(axis_stick, 1, new TGeoCombiTrans(axes_size*0.8, 0., 0., single_axis_rot)); 
        single_axis->AddNodeOverlap(axis_cone,  1, new TGeoCombiTrans(axes_size*1.8, 0., 0., single_axis_rot)); 

        auto rot_xaxis = new TGeoRotation; rot_xaxis->SetAngles(0., 0., 0.); 
        auto rot_yaxis = new TGeoRotation; rot_yaxis->SetAngles(90., 0., 0.); 
        auto rot_zaxis = new TGeoRotation; rot_zaxis->SetAngles(0., 90., 90.);

        xyz_axes->AddNodeOverlap(single_axis, 1, rot_xaxis);
        xyz_axes->AddNodeOverlap(single_axis, 2, rot_yaxis); 
        xyz_axes->AddNodeOverlap(single_axis, 3, rot_zaxis); 
    }



    TGeoRotation *rot_sieve = new TGeoRotation; 
    rot_sieve->SetAngles(90., 90., 0.); 


    //now, we're ready to add the sieve-holes
    top->AddNodeOverlap(sieve_vol, 1); 
    
    top->AddNodeOverlap(xyz_axes, 1); 

    geom->CloseGeometry(); 

    top->Draw(); 

    return 0; 

#if 0
    //--- define the transformations
    TGeoTranslation *tr1 = new TGeoTranslation(20., 0, 0.);
    TGeoTranslation *tr2 = new TGeoTranslation(10., 0., 0.);
    TGeoTranslation *tr3 = new TGeoTranslation(10., 20., 0.);
    TGeoTranslation *tr4 = new TGeoTranslation(5., 10., 0.);
    TGeoTranslation *tr5 = new TGeoTranslation(20., 0., 0.);
    TGeoTranslation *tr6 = new TGeoTranslation(-5., 0., 0.);
    TGeoTranslation *tr7 = new TGeoTranslation(7.5, 7.5, 0.);
    TGeoRotation *rot1 = new TGeoRotation("rot1", 90., 0., 90., 270., 0., 0.);
    TGeoCombiTrans *combi1 = new TGeoCombiTrans(7.5, -7.5, 0., rot1);
    TGeoTranslation *tr8 = new TGeoTranslation(7.5, -5., 0.);
    TGeoTranslation *tr9 = new TGeoTranslation(7.5, 20., 0.);
    TGeoTranslation *tr10 = new TGeoTranslation(85., 0., 0.);
    TGeoTranslation *tr11 = new TGeoTranslation(35., 0., 0.);
    TGeoTranslation *tr12 = new TGeoTranslation(-15., 0., 0.);
    TGeoTranslation *tr13 = new TGeoTranslation(-65., 0., 0.);
    
    TGeoTranslation *tr14 = new TGeoTranslation(0, 0, -100);
    TGeoCombiTrans *combi2 = new TGeoCombiTrans(0, 0, 100, new TGeoRotation("rot2", 90, 180, 90, 90, 180, 0));
    TGeoCombiTrans *combi3 = new TGeoCombiTrans(100, 0, 0, new TGeoRotation("rot3", 90, 270, 0, 0, 90, 180));
    TGeoCombiTrans *combi4 = new TGeoCombiTrans(-100, 0, 0, new TGeoRotation("rot4", 90, 90, 0, 0, 90, 0));
    TGeoCombiTrans *combi5 = new TGeoCombiTrans(0, 100, 0, new TGeoRotation("rot5", 0, 0, 90, 180, 90, 270));
    TGeoCombiTrans *combi6 = new TGeoCombiTrans(0, -100, 0, new TGeoRotation("rot6", 180, 0, 90, 180, 90, 90));
    
    //--- make the top container volume
    Double_t worldx = 110.;
    Double_t worldy = 50.;
    Double_t worldz = 5.;
    TGeoVolume *top = geom->MakeBox("TOP", Vacuum, 270., 270., 120.);
    geom->SetTopVolume(top);
    TGeoVolume *replica = geom->MakeBox("REPLICA", Vacuum, 120, 120, 120);
    replica->SetVisibility(kFALSE);
    TGeoVolume *rootbox = geom->MakeBox("ROOT", Vacuum, 110., 50., 5.);
    rootbox->SetVisibility(kFALSE);
    
    //--- make letter 'R'
    TGeoVolume *R = geom->MakeBox("R", Vacuum, 25., 25., 5.);
    R->SetVisibility(kFALSE);
    TGeoVolume *bar1 = geom->MakeBox("bar1", Al, 5., 25, 5.);
    bar1->SetLineColor(kRed);
    R->AddNode(bar1, 1, tr1);
    TGeoVolume *bar2 = geom->MakeBox("bar2", Al, 5., 5., 5.);
    bar2->SetLineColor(kRed);
    R->AddNode(bar2, 1, tr2);
    R->AddNode(bar2, 2, tr3);
    TGeoVolume *tub1 = geom->MakeTubs("tub1", Al, 5., 15., 5., 90., 270.);
    tub1->SetLineColor(kRed);
    R->AddNode(tub1, 1, tr4);
    TGeoVolume *bar3 = geom->MakeArb8("bar3", Al, 5.);
    bar3->SetLineColor(kRed);
    TGeoArb8 *arb = (TGeoArb8 *)bar3->GetShape();
    arb->SetVertex(0, 15., -5.);
    arb->SetVertex(1, 0., -25.);
    arb->SetVertex(2, -10., -25.);
    arb->SetVertex(3, 5., -5.);
    arb->SetVertex(4, 15., -5.);
    arb->SetVertex(5, 0., -25.);
    arb->SetVertex(6, -10., -25.);
    arb->SetVertex(7, 5., -5.);
    R->AddNode(bar3, 1, gGeoIdentity);
    
    //--- make letter 'O'
    TGeoVolume *O = geom->MakeBox("O", Vacuum, 25., 25., 5.);
    O->SetVisibility(kFALSE);
    TGeoVolume *bar4 = geom->MakeBox("bar4", Al, 5., 7.5, 5.);
    bar4->SetLineColor(kYellow);
    O->AddNode(bar4, 1, tr5);
    O->AddNode(bar4, 2, tr6);
    TGeoVolume *tub2 = geom->MakeTubs("tub1", Al, 7.5, 17.5, 5., 0., 180.);
    tub2->SetLineColor(kYellow);
    O->AddNode(tub2, 1, tr7);
    O->AddNode(tub2, 2, combi1);
    
    //--- make letter 'T'
    TGeoVolume *T = geom->MakeBox("T", Vacuum, 25., 25., 5.);
    T->SetVisibility(kFALSE);
    TGeoVolume *bar5 = geom->MakeBox("bar5", Al, 5., 20., 5.);
    bar5->SetLineColor(kBlue);
    T->AddNode(bar5, 1, tr8);
    TGeoVolume *bar6 = geom->MakeBox("bar6", Al, 17.5, 5., 5.);
    bar6->SetLineColor(kBlue);
    T->AddNode(bar6, 1, tr9);
    
    rootbox->AddNode(R, 1, tr10);
    rootbox->AddNode(O, 1, tr11);
    rootbox->AddNode(O, 2, tr12);
    rootbox->AddNode(T, 1, tr13);
    
    replica->AddNode(rootbox, 1, tr14);
    replica->AddNode(rootbox, 2, combi2);
    replica->AddNode(rootbox, 3, combi3);
    replica->AddNode(rootbox, 4, combi4);
    replica->AddNode(rootbox, 5, combi5);
    replica->AddNode(rootbox, 6, combi6);

    top->AddNode(replica, 1, new TGeoTranslation(-150, -150, 0));
    top->AddNode(replica, 2, new TGeoTranslation(150, -150, 0));
    top->AddNode(replica, 3, new TGeoTranslation(150, 150, 0));
    top->AddNode(replica, 4, new TGeoTranslation(-150, 150, 0));
    
    //--- close the geometry
    geom->CloseGeometry();
    
    //--- draw the ROOT box.
    // by default the picture will appear in the standard ROOT TPad.
    // if you have activated the following line in system.rootrc,
    // it will appear in the GL viewer
    // #Viewer3D.DefaultDrawOption:   ogl
    
    geom->SetVisLevel(4);
    top->Draw("");

    new TBrowser; 

    return 0; 

#endif
}