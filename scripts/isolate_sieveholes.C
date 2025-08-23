////////////////////////////////////////////////////////////////////////////////
//
//  This is a macro which creates an application which allows the creation of training data 
//  by selecting & isolating sieve-hole data 
//
//  This is the interactive app which creates a PolynomialCut, given a TH2.   
//
////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________________________________________________________
int isolate_sieveholes( const bool is_RHRS, 
                        const char* path_infile,
                        const char* coord_x ="dxdz_sv",
                        const char* coord_y ="dydz_sv",
                        const char* drawing_option ="col2", 
                        const unsigned int palette =kSunset)
{
    new PickSieveHoleApp(gClient->GetRoot(), 
                         1400, 
                         700, 
                         is_RHRS,
                         path_infile, 
                         coord_x, 
                         coord_y, 
                         drawing_option,
                         palette);

    return 0; 
}
