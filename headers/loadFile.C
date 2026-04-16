{
    ifstream vacc_file;
	// vacc_file.open("data/vacc_nojpsi.dat");
	vacc_file.open("/home/tomori/BESIII/myCodes/radInts/vaccFile/vacc_nojpsi.dat");
	 for(int i=0;i<vaccsize;i++)
		    {
			    vacc_file>>Evac[i]>>CSvac[i];
			}
    interpolator.SetData(vaccsize, Evac, CSvac);

    ifstream input_file;
    input_file.open("data/para.txt");
    int num_val=0;
    string line;
 

    while(1){
        getline(input_file,line);
            if(!input_file.good())
				break;
                // cout<<Nob[num_val]<<endl;
            input_file>>desEnergy[num_val]>>xdata[num_val]>>dEnergy[num_val]>>Nob[num_val]>>Nerr[num_val]>>epsilon[num_val]>>L[num_val]>>Lerr[num_val]>>delta[num_val]>>ddelta[num_val]>>yerrsysuncor[num_val]>>yerrsyscor[num_val];
			// input_file>>xdata[num_val]>>Nob[num_val]>>L[num_val]>>Nerr[num_val]>>epsilon[num_val]>>yerrsysuncor[num_val]>>yerrsyscor[num_val]>>delta[num_val]>>ddelta[num_val]>>dEnergy[num_val];
            num_val++;
      }
	input_file.close();
// return 0;

    // ofstream output_file;
	// TString FileName("data/profile_FF.txt");
	// output_file.open(FileName);
    for(int i=0;i<Arsize;i++)
    {
        ydata[i]= Nob[i]/epsilon[i]/L[i]/B1/B1/B2;
        yerrsta[i]= Nerr[i]/epsilon[i]/L[i]/B1/B1/B2;
        yerrsysuncor[i] = yerrsysuncor[i]*ydata[i]/100;
        yerrsyscor_quad[i] = yerrsyscor[i] * ydata[i]/100;
        yerrsys[i] = sqrt(yerrsysuncor[i]*yerrsysuncor[i]+yerrsyscor_quad[i]*yerrsyscor_quad[i]);
        yerr1[i] = sqrt(yerrsysuncor[i]*yerrsysuncor[i]+ yerrsta[i]*yerrsta[i]);
        yerr[i] = sqrt(yerrsysuncor[i]*yerrsysuncor[i]+ yerrsta[i]*yerrsta[i]+ yerrsyscor_quad[i]*yerrsyscor_quad[i]);
        // output_file<<setprecision(4)<<setiosflags(ios::fixed)<<xdata[i]<<"\t"<<ydata[i]<<"\t"<<yerrsta[i]<<"\t"<<yerrsys[i]<<"\t"<<Nob[i]<<"\t"<<epsilon[i]<<"\n";
    }
    // output_file.close();

    // vaccum polarization correction
    vpLocal.enable_vp = true;
    vaccFunc::InitVacc();
    vpLocal.vacc_func = vaccFunc::Vacc_from_table;
    vpLocal.vacc_user = vaccFunc::interp.get();

    gh128::init_once();

    std::cout << "Relative correlated systematics: " << yerrsyscor[0]/100. << std::endl;
}