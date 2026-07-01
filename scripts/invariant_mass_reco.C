#include "opticsModel/forward_production.h"
#include "include/fit_gaus_to_hist.h"
#include "include/get_n_cpus.h"
#include <PolynomialCut.h> 

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStopwatch.h> 
#include <TLegend.h> 
#include <TF1.h> 
#include <Math/LorentzVector.h>
#include <Math/Vector4D.h> 

#include <cstdio> 
#include <string> 

namespace {
    constexpr double hrs_momentum = 1104.; //MeV 
    constexpr double me2 = 0.511*0.511; //MeV^2  
}

void apply_pid_cut(RDFNodeAccumulator& rna, bool is_RHRS); 

int invariant_mass_reco(std::string path_infile, std::string path_outfile)
{   
    //get our forward-model (for production data)
    FwdProdModel model; 

    const char* tree_name = "track_data";

    using namespace ROOT::VecOps; 
    using ApexOptics::Trajectory_t; 
    using ApexOptics::SCS_to_HCS; 

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(__func__, "libApexOptics could not be loaded."); 
        return 1; 
    }

    //const bool is_RHRS = Get_TParameter_from_TFile<bool>(path_infile.data(), "is_RHRS").value(); 
    
    const unsigned int n_cpus = get_n_cpus(); 
    ROOT::EnableImplicitMT(n_cpus); 
    Info(__func__, "Multi-threadding is enabled. cpus available: %u", n_cpus); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //fit the background DT spectrum 
    auto hist_dt_fit = df
        
        .Filter([](const RVec<double>& R, const RVec<double>& L){ 
            return R.size() == L.size() == 1; 
        }, {"R_tracks_S2time", "L_tracks_S2time"})
        
        .Define("R_L_dt", [](const RVec<double>& R, const RVec<double>& L){
            return (R[0] - L[0])/1e-9; 
        }, {"R_tracks_S2time", "L_tracks_S2time"})

        .Histo1D<double>({"h_dt", ";T_{R} - T_{L} (ns);", 200, 0., 80}, "R_L_dt");

    double center, sigma, amplitude, base; 

    fit_gaus_to_hist((TH1D*)hist_dt_fit->Clone("dt_clone"), 5., center, sigma, amplitude, base, true);
    
    std::printf(
        "gauss fit results:\n"
        "   center: %-4.1f ns\n"
        "   sigma:  %-4.3f ns\n"
        "   signal max amplitude: %-7.0f\n"
        "   background amplitude: %-7.0f\n",
        center, sigma, amplitude, base
    );


    ULong64_t n_events_total = *df.Count(); 
    
    //this is a little helper class which is meant to avoid some of the awkward syntax typically associated with RDataFrames creation. 
    auto rna = RDFNodeAccumulator(df);  

    //define some parameters about the signal/noise ratio 
    rna.Define("R_L_dt", [](const RVec<double>& R, const RVec<double>& L){
            return (R[0] - L[0])/1e-9; 
        }, {"R_tracks_S2time", "L_tracks_S2time"});

    double max_snr_stddev = 5.;
    rna = rna.Get().Filter([max_snr_stddev, sigma, center](double R_L_dt){ 
        return std::fabs(R_L_dt - center)/sigma < max_snr_stddev; 
    }, {"R_L_dt"});

    //prob. that this event is a 'true' coinc. event
    rna.Define("p_coinc", [sigma, center, amplitude, base](double R_L_dt){
            double x = (R_L_dt - center)/sigma;

            //compute relative signal / bg amplitudes 
            double s = amplitude * std::exp( -x*x/2. );
            double b = base; 
            
            return s / (s + b); 
            
        }, {"R_L_dt"});

        

    //define output branches
    rna = model.DefineOutputs(rna.Get()); 
    
    //do PID cuts
    apply_pid_cut(rna, true); 
    apply_pid_cut(rna, false); 

    using FourVec = ROOT::Math::XYZTVector; 

    //define invariant mass
    //sum of 3-momenta of both positron and electron 
    rna.Define("P_p", [](const Trajectory_t& R_Xsv_scs)
    {
        auto R_Xsv = SCS_to_HCS(true,  R_Xsv_scs);
        TVector3 R_p( R_Xsv.dxdz, R_Xsv.dydz, 1. ); R_p = R_p.Unit() * hrs_momentum*( 1. + R_Xsv.dpp );
        return R_p; 
    }, {"R_Xsv_reco"}); 

    rna.Define("P_e", [](const Trajectory_t& L_Xsv_scs)
    {
        auto L_Xsv = SCS_to_HCS(false, L_Xsv_scs);
        TVector3 L_p( L_Xsv.dxdz, L_Xsv.dydz, 1. ); L_p = L_p.Unit() * hrs_momentum*( 1. + L_Xsv.dpp );
        return L_p; 
    }, {"L_Xsv_reco"}); 

    rna.Define("P_sum", [](const TVector3& Pe, const TVector3& Pp)
    {
        auto sum = Pp + Pe; 

        return FourVec( 
            sum.x(), 
            sum.y(), 
            sum.z(),
            std::sqrt( Pp.Mag2() + me2 ) + std::sqrt( Pe.Mag2() + me2 )
        ); 
        
    }, {"P_e", "P_p"}); 

    rna.Define("reco_invariant_mass", [](const FourVec& P){ return std::sqrt( P.M2() ); }, {"P_sum"}); 

    rna.Define("x_hcs", [](TVector3 vtx){return vtx.x();}, {"position_vtx_reco"});
    rna.Define("z_hcs", [](TVector3 vtx){return vtx.z();}, {"position_vtx_reco"});

    rna = rna.Get().Filter([](double x_hcs, double z_hcs){
        if (z_hcs > 0.300 || z_hcs < -0.300) return false; 
        if (-0.02 > x_hcs || x_hcs > 0.04) return false; 
        return true; 
    }, {"x_hcs", "z_hcs"}); 

    std::vector<std::string> R_cut_path_list{
        "data/polycuts/RHRS_dx_dp.dat",
        "data/polycuts/RHRS_dy_dp.dat",
        "data/polycuts/RHRS_dx_dy.dat"
    }; 
    std::vector<std::string> L_cut_path_list{
        "data/polycuts/LHRS_dx_dp.dat",
        "data/polycuts/LHRS_dy_dp.dat",
        "data/polycuts/LHRS_dx_dy.dat"
    }; 
    std::vector<PolynomialCut> polycuts_R; polycuts_R.reserve(R_cut_path_list.size()); 
    std::vector<PolynomialCut> polycuts_L; polycuts_L.reserve(L_cut_path_list.size()); 
    
    for (auto& cut_path : R_cut_path_list) {
        polycuts_R.emplace_back(); 
        auto& polycut = polycuts_R.back(); 
        polycut.Parse_dbfile(cut_path.c_str()); 
    }
    for (auto& cut_path : L_cut_path_list) {
        polycuts_L.emplace_back(); 
        auto& polycut = polycuts_L.back(); 
        polycut.Parse_dbfile(cut_path.c_str()); 
    }
    rna = rna.Get()

        //filter right arm     
        .Filter([&polycuts_R](const Trajectory_t& t){ 
            if (polycuts_R[0].IsInside(t.dxdz, t.dpp) &&
                polycuts_R[1].IsInside(t.dydz, t.dpp) &&
                polycuts_R[2].IsInside(t.dxdz, t.dydz)) return true; 
            
            return false; 
        }, {"R_Xsv_reco"})

        //filter left arm 
        .Filter([&polycuts_L](const Trajectory_t& t){ 
            if (polycuts_L[0].IsInside(t.dxdz, t.dpp) &&
                polycuts_L[1].IsInside(t.dydz, t.dpp) &&
                polycuts_L[2].IsInside(t.dxdz, t.dydz)) return true; 
            
            return false; 
        }, {"L_Xsv_reco"});


    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(__func__, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    //create both histograms
    auto R_hist_xy = rna.Get()
        .Histo2D({"h_xy",     "(RHRS) Sieve-plane projection;x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, +0.025}, "R_x_sv", "R_y_sv");
    
    auto R_hist_angles = rna.Get()
        .Histo2D({"h_R_angles", "(RHRS) Sieve-plane projection;dx/dx_{sv};dy/dz_{sv}", 200, -0.065,0.065, 200, -0.065,+0.065}, "R_dxdz_sv", "R_dydz_sv"); 
    
    auto R_hist_dxdz_dp = rna.Get()
        .Histo2D({"h_R_dx_dp", "(RHRS) Sieve-plane projection;dx/dx_{sv};dp/p_{sv}", 200, -0.065,+0.065, 200, -0.06,+0.06}, "R_dxdz_sv", "R_dpp_sv"); 
    
    auto R_hist_dydz_dp = rna.Get()
        .Histo2D({"h_R_dy_dp", "(RHRS) Sieve-plane projection;dy/dx_{sv};dp/p_{sv}", 200, -0.06,+0.06, 200, -0.06,+0.06}, "R_dydz_sv", "R_dpp_sv"); 
    

    //create both histograms
    auto L_hist_xy = rna.Get()
        .Histo2D({"h_xy",     "(LHRS) Sieve-plane projection;x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, +0.025}, "L_x_sv", "L_y_sv");
    
    auto L_hist_angles = rna.Get()
        .Histo2D({"h_L_angles", "(LHRS) Sieve-plane projection;dx/dx_{sv};dy/dz_{sv}", 200, -0.065,+0.065, 200, -0.065,+0.065}, "L_dxdz_sv", "L_dydz_sv"); 
    
    auto L_hist_dxdz_dp = rna.Get()
        .Histo2D({"h_L_dx_dp", "(LHRS) Sieve-plane projection;dx/dx_{sv};dp/p_{sv}", 200, -0.065,+0.065, 200, -0.06,+0.06}, "L_dxdz_sv", "L_dpp_sv"); 
    
    auto L_hist_dydz_dp = rna.Get()
        .Histo2D({"h_L_dy_dp", "(LHRS) Sieve-plane projection;dy/dx_{sv};dp/p_{sv}", 200, -0.06,+0.06, 200, -0.06,+0.06}, "L_dydz_sv", "L_dpp_sv"); 
    
    auto his_xz_hcs = rna.Get()
        .Histo2D({"h_xz", "React vertex reconstruction (x & z);z_{hcs} (m);x_{hcs} (m)", 200, -0.30,+0.30, 200,-0.30,+0.30}, "z_hcs", "x_hcs"); 

    auto hist_m = rna.Get() 
        .Histo1D<double>({"h_m", "Invariant mass (MeV);m_{A} (MeV);counts / MeV", 200, 100, 300}, "reco_invariant_mass"); 

    auto hist_dt = rna.Get()
        .Histo1D<double>({"h_dt", "T_{R} - T_{L} (ns)", 
            45, center - sigma*max_snr_stddev, center + sigma*max_snr_stddev}, "R_L_dt"); 

    auto hist_m_accident = rna.Get()
        .Filter([](double pc){ return pc < 0.05; }, {"p_coinc"})
        .Histo1D<double>({"h_m", "Invariant mass (MeV);m_{A} (MeV);counts / MeV", 180, 120, 300}, "reco_invariant_mass"); 
    
    auto hist_m_signal = rna.Get()
        .Filter([](double pc){ return pc > 0.8; }, {"p_coinc"})
        .Histo1D<double>({"h_m", "Invariant mass (MeV);m_{A} (MeV);counts / MeV", 180, 120, 300}, "reco_invariant_mass"); 

    //sum of 3-vector of both positron and electron 
    rna.Get().Snapshot("track_data", path_outfile, {
        "P_e", 
        "P_p", 
        "position_vtx_reco", 
        "vertex_closest_approach_distance", 
        "p_coinc"
    });

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    //gStyle->SetOptStat(0); 

    //PolynomialCut::InteractiveApp((TH2*)hist_z_y->Clone("hclone"), "col", kSunset); 
    //return 0; 
    
    auto outfile = new TFile(path_outfile.c_str(), "UPDATE"); 

    TCanvas* c; 
    
    c = new TCanvas("c1_R", Form("data: %s",path_infile.data()), 1200, 600); 
    gPad->SetLeftMargin(0.15); 
    c->Divide(2,1, 0.01,0.01); 

    TStopwatch timer; 

    (c->cd(1))->SetLeftMargin(0.15); { R_hist_angles->DrawCopy("col"); R_hist_angles->Write(); }
    (c->cd(2))->SetLeftMargin(0.15); { L_hist_angles->DrawCopy("col"); L_hist_angles->Write(); }

    double user_time = timer.CpuTime(); 
    std::printf(" total user time: %.2f s. ( %.3f us / event )\n", user_time, 1e6*user_time/((double)n_events_total));


    c = new TCanvas("c1_L_angles_dp", Form("data: %s",path_infile.data()), 1200, 600); 
    gPad->SetLeftMargin(0.15); 
    c->Divide(2,1, 0.01,0.01); 

    (c->cd(1))->SetLeftMargin(0.15); { L_hist_dxdz_dp->DrawCopy("col"); L_hist_dxdz_dp->Write(); } 
    (c->cd(2))->SetLeftMargin(0.15); { L_hist_dydz_dp->DrawCopy("col"); L_hist_dydz_dp->Write(); }


    c = new TCanvas("c1_R_angles_dp", Form("data: %s",path_infile.data()), 1200, 600); 
    gPad->SetLeftMargin(0.15); 
    c->Divide(2,1, 0.01,0.01); 

    (c->cd(1))->SetLeftMargin(0.15); { R_hist_dxdz_dp->DrawCopy("col"); R_hist_dxdz_dp->Write(); } 
    (c->cd(2))->SetLeftMargin(0.15); { R_hist_dydz_dp->DrawCopy("col"); R_hist_dydz_dp->Write(); }


    c = new TCanvas("c1_xz", Form("data: %s",path_infile.data()), 1200, 600); 
    gPad->SetLeftMargin(0.15); 

    his_xz_hcs->DrawCopy("col"); 
    his_xz_hcs->Write(); 

    c = new TCanvas("c_m", Form("data: %s",path_infile.data()), 800, 600);
    hist_m->DrawCopy(); 
    hist_m->Write(); 


    c = new TCanvas; 
    hist_dt->SetMaximum( hist_dt->GetMaximum()*1.15 );
    hist_dt->SetMinimum( 0. );
    hist_dt->DrawCopy("HIST, E"); 
    hist_dt->Write(); 

    fit_gaus_to_hist((TH1D*)hist_dt->Clone("dt_clone2"), 4., center, sigma, amplitude, base, false);
    
    std::printf(
        "new, post-filtering gauss fit results:\n"
        "   center: %-4.1f ns\n"
        "   sigma:  %-4.3f ns\n"
        "   signal max amplitude: %-7.0f\n"
        "   background amplitude: %-7.0f\n",
        center, sigma, amplitude, base
    );

    auto tf1 = new TF1("fcn", [center, sigma, amplitude, base](double *x, double *par){
        double arg = (x[0] - center)/sigma; 
        return base + amplitude*std::exp( -arg*arg/2. ); 
    },
        hist_dt->GetXaxis()->GetXmin(), 
        hist_dt->GetXaxis()->GetXmax(), 
        0
    ); 

    tf1->Draw("SAME"); 

    //prob. that this event is an accidental coinc. event
    auto hist_m_p = rna.Get()
    
        .Redefine("p_coinc", [sigma, center, amplitude, base](double R_L_dt){
            double x = (R_L_dt - center)/sigma;

            //compute relative signal / bg amplitudes 
            double s = amplitude * std::exp( -x*x/2. );
            double b = base; 
            
            return s / (s + b); 
            
        }, {"R_L_dt"})

        .Histo2D<double>({"h_m_p", ";m_{#pm} (MeV);p. accidental", 180, 120, 300, 50, 0., 1.}, "reco_invariant_mass", "p_coinc"); 
    
    c = new TCanvas; 
    hist_m_p->DrawCopy("col");
    hist_m_p->Write(); 

    c = new TCanvas; 
    auto legend = new TLegend; 

    auto hm_a = (TH1D*)hist_m_accident->Clone("accident");
    hm_a->SetLineColor(kRed); 
    hm_a->SetFillColor(kRed);
    hm_a->SetFillStyle(3004); 
    hm_a->Scale( 1./hm_a->GetMaximum() );
    hm_a->Draw("HIST"); 
    legend->AddEntry(hm_a, "accidental"); 

    auto hm_s = (TH1D*)hist_m_signal->Clone("signal");
    hm_s->SetLineColor(kBlack);
    hm_s->Scale( 1./hm_s->GetMaximum() );
    hm_s->Draw("SAME HIST"); 
    legend->AddEntry(hm_s, "signal"); 

    legend->Draw(); 

    outfile->Close(); 

    return 0; 
}

void apply_pid_cut(RDFNodeAccumulator& rna, bool is_RHRS)
{
    using ApexOptics::Trajectory_t; 

    //minimum sum of all cerenkov ADCs. 
    const double min_cerenkov_sum = is_RHRS ? 1.0e3 : 1.5e3;  
    const double min_Esh_Eps_sum  = is_RHRS ? 0.45  : 0.45; 

    std::string branch_cerenkov_adc     = is_RHRS ? "R_cer_adc" : "L_cer_adc"; 
    std::string branch_preshower_adc    = is_RHRS ? "R_ps_adc" : "L_ps_adc"; 
    std::string branch_shower_adc       = is_RHRS ? "R_sh_adc" : "L_sh_adc"; 
    
    std::string arm = is_RHRS ? "R" : "L"; 

    //perform a very basic cut on sum of cerenkov ADCs 
    rna.Define(arm+"_cer_sum", [](const ROOT::RVec<double>& cer_a_c)
        {
            double val=0.; 
            for (double adc : cer_a_c) if (fabs(adc) < 5e4) val += adc; 
            return val; 
        }, {branch_cerenkov_adc});

        //for now, this is actually only the MAXIMUM block. 
    rna.Define(arm+"_E_ps", [](const ROOT::RVec<double>& blocks_adc)
        {
            double val=0.; 
            for (double adc : blocks_adc) val = std::max(val, adc);
            return val; 
        }, {branch_preshower_adc});

        //for now, this is actually only the MAXIMUM block. 
    rna.Define(arm+"_E_sh", [](const ROOT::RVec<double>& blocks_adc)
        {
            double val=0.; 
            for (double adc : blocks_adc) val = std::max(val, adc); //val += adc; 
            return val; 
        }, {branch_shower_adc}); 

       
        //ratio of shower energy over momentum 
    rna.Define(arm+"_E_sh_p_ratio", [](double E, Trajectory_t Xsv)
        {
            return E / (hrs_momentum * ( 1. + Xsv.dpp ));  
        }, {arm+"_E_sh", arm+"_Xsv_reco"}); 

        //ratio of pre-shower energy over momentum 
    rna.Define(arm+"_E_ps_p_ratio", [](double E, Trajectory_t Xsv)
        {
            return E / (hrs_momentum * ( 1. + Xsv.dpp ));  
        }, {arm+"_E_ps", arm+"_Xsv_reco"}); 

    rna = rna.Get().Filter([min_Esh_Eps_sum, min_cerenkov_sum](double Eps, double Esh, double cer_sum){
        return (Eps + Esh > min_Esh_Eps_sum) && (cer_sum > min_cerenkov_sum); 
    }, {arm+"_E_ps_p_ratio", arm+"_E_sh_p_ratio", arm+"_cer_sum"}); 

    return; 
}