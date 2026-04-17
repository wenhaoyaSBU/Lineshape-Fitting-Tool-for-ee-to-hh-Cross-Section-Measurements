#include "physicsFuncs.h"
#include "variables.h"
#include "BrErr.h"

int main()
{
     
	auto start = std::chrono::high_resolution_clock::now();
 
    #include "loadFile.C"  // Load data from file and do some pre-calculation


    // save the calculated cross section data (ydata) to a file
    ofstream xs_data_file;
    TString xsFileName("output/calculated_xs_data.txt");
    xs_data_file.open(xsFileName);
    xs_data_file<<"#Energy(GeV)\tCrossSection(nb)\tStatError(nb)\n";
    for(int i=0;i<Arsize;i++)
    {
        xs_data_file<<setprecision(6)<<setiosflags(ios::fixed)<<xdata[i]<<"\t"<<ydata[i]<<"\t"<<yerrsta[i]<<"\n";
    }
    xs_data_file.close();
    
 

	TMinuit *gMinuit = new TMinuit(Param+1);
    cout<<"fcn "<<&fcn<<endl;
    gMinuit->SetFCN(fcn);
    // cout<<"fcn "<<&fcn_linearW<<endl;
    // gMinuit->SetFCN(fcn_linearW);

    Int_t ierflg = 0;
    double arglist[10];

    
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    const int pa1 = numPara;  //非pull项的个数


    // double vstart0[pa1] = {3.09688, 1.2622, 0, 4.43e-1, -8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // negative solution
    // double vstart0[pa1] = {3.09688, 1.2622, 0, 4.43e-1, 8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // positive solution
    // double vstart0[pa1] = {3.09688, 1.2622, 0, 1.38, 8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // positive solution
    
    double vstart0[pa1] = {MJpsi_mean, 1.2622, 0, 7.8, 8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // positive solution, W^-8
    // double vstart0[pa1] = {MJpsi_mean, 1.2622, 0, 4.23, 8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // positive solution, W^-7
    // double vstart0[pa1] = {MJpsi_mean, 1.2622, 0, 2.41, 8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // positive solution, W^-6

    // double vstart0[pa1] = {MJpsi_mean, 1.2622, 0, 7.5, -8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // negative solution, W^-8
    // double vstart0[pa1] = {MJpsi_mean, 1.2622, 0, 4.23, -8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // negative solution, W^-7
    // double vstart0[pa1] = {MJpsi_mean, 1.2622, 0, 2.41, -8.35575e-1, 0, Gamma_mean, Gamma_ee_mean, 0.00093, 1.};       // negative solution, W^-6
    
    
    // double vstart0[pa1] = {3.0969, 1.15, 0, 4.8e-1, PI/2, 0, 0.0008, 1.};
    // double vstart0[pa1] = {3.097, 1.15, 0, 4.8e-1, 0, 0, 0.0008, 1.};
    // double vstart0[pa1] = {3.097, 1.00, 0, 4.8e-1, 0, 0, 0.0008, 1.};

    double step0[pa1]   = { 1e-3,1e-3,1e-3,1e-5, 1e-3, 1e-3, Gamma_std*0.1, Gamma_ee_std*0.1, 1e-7, 1e-3};
 
    // double vstart0[pa1] = {3.097, 2.50, 0, 2.12e-1, 0, 0, 0.0008};

    // double step0[pa1]   = { 1e-3,1e-3,1e-3,1e-5, 1e-3, 1e-3, 1e-7};
 
 //   double step0[4]   = { 0.0,0.0,0.0,0.0};

    for(int i=0;i<pa1;i++)
    {
        vstart[i]=vstart0[i];
        step[i]=step0[i];
 //       step[i]=0.0;
    }

    for(int i=0;i<Arsize;i++)
    {
        vstart[i+pa1]=xdata[i];
 //       step[i+4]=1e-9;
        step[i+pa1]=dEnergy[i]/100.;
    }

    gMinuit->mnparm(0, "M", vstart[0], step[0], 3.09,    3.10,   ierflg);	
    gMinuit->mnparm(1, "CC1", vstart[1], step[1], 0.8,    2.2,   ierflg);
    gMinuit->mnparm(2, "CC2 ", vstart[2], step[2], 0,    60,   ierflg);
    // gMinuit->mnparm(3, "FF ", vstart[3], step[3], 0.3,     0.55,   ierflg);      // W^-6
    // gMinuit->mnparm(3, "FF ", vstart[3], step[3], 0.0,     4.5,   ierflg);
    gMinuit->mnparm(3, "FF ", vstart[3], step[3], 4.8,     10.0,   ierflg);       // W^-8
    // gMinuit->mnparm(3, "FF ", vstart[3], step[3], 2.74,     5.71,   ierflg);       // W^-7
    // gMinuit->mnparm(3, "FF ", vstart[3], step[3], 1.5,     3.3,   ierflg);       // W^-6
    gMinuit->mnparm(4, "phi1", vstart[4], step[4], 0,    PI,   ierflg);	//positive
    // gMinuit->mnparm(4, "phi1", vstart[4], step[4], -PI,    0,   ierflg);	//negative
    // gMinuit->mnparm(4, "phi1", vstart[4], step[4], -PI,    PI,   ierflg);	//all
	gMinuit->mnparm(5, "phi2", vstart[5], step[5], -PI,    PI,   ierflg);
	gMinuit->mnparm(6, "Gamma", vstart[6], step[6], Gamma_mean - 10*Gamma_std,    Gamma_mean + 10*Gamma_std,   ierflg);
	gMinuit->mnparm(7, "Gamma_ee", vstart[7], step[7], Gamma_ee_mean - 10*Gamma_ee_std,    Gamma_ee_mean + 10*Gamma_ee_std,   ierflg);
    gMinuit->mnparm(8, "SE", vstart[8], step[8], 0.00079,    0.00115,   ierflg);
    gMinuit->mnparm(9, "f", vstart[9], step[9], 0.8,    1.2,   ierflg);

	gMinuit->FixParameter(2);
	// gMinuit->FixParameter(1);
    gMinuit->FixParameter(5);

    // Fix M, Gamma and Gamma_ee
    // gMinuit->FixParameter(0);
    // gMinuit->FixParameter(6);
    // gMinuit->FixParameter(7);

    // PURE EM
    // gMinuit->FixParameter(1);
    // gMinuit->FixParameter(4);

    if(yerrsyscor[0] == 0) gMinuit->FixParameter(7);

    string nameM[Arsize];

    for(int i=0;i<Arsize;i++)
    {
        nameM[i]="M"+to_string(i);
        gMinuit->mnparm(i+pa1, nameM[i], vstart[i+pa1], step[i+pa1], xdata[i]-15*dEnergy[i],    xdata[i]+15*dEnergy[i],   ierflg);

        // For linear approx.
        // gMinuit->FixParameter(i+pa1);
    }

    arglist[0] = 1000000;
    arglist[1] = 0.001;
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    gMinuit->mnexcm("HESSE", arglist, 2, ierflg);

    double amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    gMinuit->mnprin(0,amin);


    double fittedParr[numPara];
    double fittedParrErr[numPara];
    for(int i = 0; i < numPara; i++){
        double paraValue, paraError;
        gMinuit->GetParameter(i, paraValue, paraError);
        fittedParr[i] = paraValue;
        fittedParrErr[i] = paraError;
    }
    ofstream outputFile("output/getpoint.txt");
    double startW = 3.05;
    double endW = 3.12;
    int nW = 2000;
    double step = (endW - startW) / (nW-1);
    for(double i=startW;i<endW;i=i+step)
    {
         outputFile <<setprecision(15)<< i <<"\t"<<Conv(i,fittedParr)<<endl;
    }
    // other components
    ofstream outputFile_con("output/getpoint_con.txt");
    for(double i=startW;i<endW;i=i+step)
    {
         outputFile_con <<setprecision(15)<< i <<"\t"<<Conv_component(i,fittedParr, "con")<<endl;
    }
    ofstream outputFile_res("output/getpoint_res.txt");
    for(double i=startW;i<endW;i=i+step)
    {
         outputFile_res <<setprecision(15)<< i <<"\t"<<Conv_component(i,fittedParr, "res")<<endl;
    }
    ofstream outputFile_int("output/getpoint_int.txt");
    for(double i=startW;i<endW;i=i+step)
    {
         outputFile_int <<setprecision(15)<< i <<"\t"<<Conv_component(i,fittedParr, "int")<<endl;
    }
    

    // points txt file output
    ofstream dressedXSFile("output/getpoint_dressed.txt");
    for(double i=startW;i<endW;i=i+step)
    {
        double dressedValue = sigma_born_dressed(i, fittedParr[4], fittedParr[3], fittedParr[1], interpolator.Eval(i), fittedParr[0], "total", true);
        dressedXSFile <<setprecision(15)<< i <<"\t"<<dressedValue<<endl;
    }
    ofstream dressedXSFile_con("output/getpoint_dressed_con.txt");
    for(double i=startW;i<endW;i=i+step)
    {
        double dressedValue = sigma_born_dressed(i, fittedParr[4], fittedParr[3], fittedParr[1], interpolator.Eval(i), fittedParr[0], "con", true);
        dressedXSFile_con <<setprecision(15)<< i <<"\t"<<dressedValue<<endl;
    }
    ofstream dressedXSFile_res("output/getpoint_dressed_res.txt");
    for(double i=startW;i<endW;i=i+step)
    {
        double dressedValue = sigma_born_dressed(i, fittedParr[4], fittedParr[3], fittedParr[1], interpolator.Eval(i), fittedParr[0], "res", true);
        dressedXSFile_res <<setprecision(15)<< i <<"\t"<<dressedValue<<endl;
    }
    ofstream dressedXSFile_int("output/getpoint_dressed_int.txt");
    for(double i=startW;i<endW;i=i+step)
    {
        double dressedValue = sigma_born_dressed(i, fittedParr[4], fittedParr[3], fittedParr[1], interpolator.Eval(i), fittedParr[0], "int", true);
        dressedXSFile_int <<setprecision(15)<< i <<"\t"<<dressedValue<<endl;
    }
    ofstream dressedXSFile_intalt("output/getpoint_dressed_int_alt.txt");
    for(double i=startW;i<endW;i=i+step)
    {
        double dressedValue = sigma_born_dressed(i, fittedParr[4], fittedParr[3], fittedParr[1], interpolator.Eval(i), fittedParr[0], "int_alt", true);
        dressedXSFile_intalt <<setprecision(15)<< i <<"\t"<<dressedValue<<endl;
    }
    // gammagg component
    if(fittedParr[1] > 0){
        double Cgammagg = fittedParr[1] - 1;
        ofstream dressedXSFile_res_gammagg("output/getpoint_dressed_res_gammagg.txt");
        for(double i=startW;i<endW;i=i+step)
        {
            double dressedValue = sigma_born_dressed(i, fittedParr[4], fittedParr[3], Cgammagg, interpolator.Eval(i), fittedParr[0], "res", true);
            dressedXSFile_res_gammagg <<setprecision(15)<< i <<"\t"<<dressedValue<<endl;
        }
    }


     


    ofstream outFilePara("output/fittedParameters.txt");
    for(int i = 0; i < numPara; i++){
        outFilePara << setprecision(15) << fittedParr[i] << "\t" << fittedParrErr[i] << endl;
    }
    // Print the fitted chi^2 here
    outFilePara << setprecision(15) << amin << "\t" << 0.0 << endl;

    ofstream outFileM("output/fittedM.txt");
    for(int i = 0; i < Arsize; i++){
        double fittedMvalue, fittedMerror;
        double fittedMshift;
        gMinuit->GetParameter(numPara + i, fittedMvalue, fittedMerror);
        fittedMshift = fittedMvalue - xdata[i];
        outFileM << setprecision(15) << fittedMvalue << "\t" << fittedMerror << "\t" << desEnergy[i] << "\t" << dEnergy[i] << "\t" << fittedMshift << endl;
    }

    // Calculate Br and BrEM with propagated uncertainties
    BrErr::Result brRes   = BrErr::CalcBrAndErr(gMinuit, fittedParr, numPara);
    BrErr::Result brEMRes = BrErr::CalcBrEMAndErr(gMinuit, fittedParr, numPara);
    
    double Brpost    = brRes.value;
    double Brerr     = brRes.error;
    double BrEMpost  = brEMRes.value;
    double BrEMerr   = brEMRes.error;
    
    std::cout << std::setprecision(15)
              << "Br   = " << Brpost   << " +- " << Brerr   << std::endl;
    std::cout << std::setprecision(15)
              << "BrEM = " << BrEMpost << " +- " << BrEMerr << std::endl;

    // Open a file and save the Br and BrEM results
    ofstream brFile("output/fitted_Br.txt");
    brFile << std::setprecision(15)
            << Brpost   << " \t " << Brerr   << std::endl;
    brFile << std::setprecision(15)
            << BrEMpost << " \t " << BrEMerr << std::endl;
    

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Running time: " << ": " << duration.count() / 1000000. << " seconds" << std::endl;

    return 0;
}
