#include <iostream>
#include <iomanip>
#include <fstream>
#include "phydro.h"
using namespace std;


double err(double x, double ref){
	return std::min(abs(x-ref), abs(x/ref-1));
}

int check(double x, double ref, double err=1e-6){
	//cout << "err: " << abs(x-ref) << " " << abs(x/ref-1)<< "\n";
	//cout << "comp: " << x << " " << ref << "|" << abs(x-ref) << "\n";
	if (abs(x-ref) < err || abs(x/ref-1) < err) return 0;
	else return 1;
}


vector<double> seq(double start, double end, int length){
	vector<double> x;
	for (int i=0; i<length; ++i) x.push_back(start + double(i)*(end-start)/(length-1));
	return x;
}

vector<double> lseq(double start, double end, int length){
	vector<double> x;
	for (int i=0; i<length; ++i) x.push_back(exp(log(start) + double(i)*(log(end)-log(start))/(length-1)));
	return x;
}

int main(){
	vector<double> x = seq(-5,40,46);
	// ftemp_inst_vcmax
	vector<double> expected_kthunge = {0.16224749762063, 0.173636341666431, 0.185731247018254, 0.19856989574356, 0.212191620359412, 0.226637452541204, 0.241950169605546, 0.258174337337494, 0.275356347157117, 0.293544444822553, 0.312788746763373, 0.333141238616272, 0.354655748442466, 0.377387884237216, 0.401394921421314, 0.426735620665162, 0.453469949152501, 0.481658668605213, 0.511362740234979, 0.542642479202111, 0.575556367816386, 0.610159405985626, 0.646500837405835, 0.684621038620159, 0.724547293377837, 0.766288095388894, 0.809825529038732, 0.855105174035919, 0.902022877428985, 0.950407658283492, 1, 1.05042491711027, 1.10115957532813, 1.15149606748754, 1.20050142672755, 1.2469793282723, 1.28944133406019, 1.32609979559137, 1.35489880528583, 1.37360195980551, 1.37995303612314, 1.37191431817576, 1.34796519490302, 1.30741376997327, 1.25064729378096, 1.17923975686556};

	for (size_t i=0; i<x.size(); ++i){
		double T = x[i];
		double f = phydro::calc_ftemp_inst_jmax(T, 25, 25, 25, phydro::FV_kumarathunge19);
		if (check(f, expected_kthunge[i]) == 1) return 1;
		//else cout << f << " " << expected_kthunge[i] << "\n";
	}


	vector<double> expected_kattge = {0.108165296210614, 0.1175411330303, 0.127651024795746, 0.138546040580533, 0.150280321556271, 0.162911228097918, 0.176499488675189, 0.19110934874336, 0.20680871702301, 0.2236693053978, 0.241766757029895, 0.261180755011782, 0.281995100691034, 0.30429774637993, 0.328180761031501, 0.353740199003026, 0.381075830401281, 0.410290675603626, 0.441490264927, 0.474781515212076, 0.510271075968114, 0.548062945862371, 0.588255092510554, 0.630934721404157, 0.676171730670059, 0.724009756496315, 0.774454063444923, 0.82745537775353, 0.882888630174217, 0.940525526345156, 1, 1.06076609384524, 1.12204890000608, 1.18279119576327, 1.24160166415102, 1.29671527092002, 1.3459821620838, 1.38690701584358, 1.41676324547863, 1.43280134465593, 1.4325533122752, 1.41420402895972, 1.37696211598528, 1.32133373738511, 1.24920421515943, 1.16367598004045};

	for (size_t i=0; i<x.size(); ++i){
		double T = x[i];
		double f = phydro::calc_ftemp_inst_jmax(T, 25, 25, 25, phydro::FV_kattge07);
		if (check(f, expected_kattge[i]) == 1) return 1;
		//else cout << f << " " << expected_kattge[i] << "\n";
	}


	vector<double> expected_leuning = {0.120375999968243, 0.130894617257924, 0.142242095547883, 0.154475762525289, 0.167656087918668, 0.181846711992658, 0.197114428717115, 0.213529106725458, 0.231163525951592, 0.250093101138906, 0.270395454897424, 0.292149792245117, 0.315436015172373, 0.340333499283094, 0.366919434633578, 0.395266609349475, 0.425440487691548, 0.457495404927633, 0.491469671831667, 0.527379355975122, 0.565210492197857, 0.604909481860986, 0.646371486311908, 0.68942672777097, 0.733824811005138, 0.779217507449454, 0.825140935113514, 0.870998744771333, 0.916048772557816, 0.959396560997767, 1, 1.03668977973091, 1.06820994072775, 1.09328108237436, 1.11068543953201, 1.11936819194594, 1.1185438676821, 1.10779211771609, 1.08712541568215, 1.05701387456219, 1.01835951405898, 0.972422293501319, 0.920710029094153, 0.864850840297628, 0.806468217168184, 0.747075353811709};

	for (size_t i=0; i<x.size(); ++i){
		double T = x[i];
		double f = phydro::calc_ftemp_inst_jmax(T, 25, 25, 25, phydro::FV_leuning02);
		if (check(f, expected_leuning[i]) == 1) return 1;
		//else cout << f << " " << expected_leuning[i] << "\n";
	}

	return 0;

}
