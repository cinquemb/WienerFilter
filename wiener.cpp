#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>


std::vector<double> cross_correlation_wiener(std::vector<double>& ts, int& new_order, int& old_order){
	//see: http://docs.scipy.org/doc/numpy/reference/generated/numpy.correlate.html
	std::vector<double> out_data;
	int n_taps = ts.size();
	int n_size;
	if(n_taps > new_order)
		n_size = n_taps;
	else
		n_size = new_order;

	for(int i = -n_taps; i< n_taps; ++i){
		double ts_sum = 0.0;
		for(int j=0; j< n_size; ++j){
			if( ((j+i) < n_taps) && (j < new_order) && ((j+i) > -1) )
				ts_sum += ts[j+i];
		}
		out_data.push_back(ts_sum);
	}

	std::vector<double> out(out_data.end()-old_order-n_taps, out_data.end()-old_order);
	return out;
}

std::vector<double> wiener_filter_1d(std::vector<double>& t_series, int& order){
	//see def wiener: https://github.com/scipy/scipy/blob/master/scipy/signal/signaltools.py
	int new_order = (2*order)+1;
	std::vector<double> t_2_series, lmean, lvar, res;

	for(auto it = t_series.begin(); it != t_series.end(); ++it)
		t_2_series.push_back(std::pow(*it, 2));

	//Estimate the local mean
	std::vector<double> t_lmean = cross_correlation_wiener(t_series, new_order, order);
	for(auto it = t_lmean.begin(); it != t_lmean.end(); ++it)
		lmean.push_back((*it)/(double)new_order);
	
	//Estimate the local variance
	std::vector<double> t_lvar = cross_correlation_wiener(t_2_series, new_order, order);
	for(int i=0; i<t_lvar.size(); ++i)
		lvar.push_back(((t_lvar[i])/(double)new_order) - (std::pow(lmean[i], 2)));

	//Estimate the noise power
	double est_n = (double)std::accumulate(lvar.begin(), lvar.end(), 0.0) / (double)lvar.size();

	for(int i=0; i<t_series.size(); ++i){
		if (lvar[i] < est_n)
			res.push_back(lmean[i]);
		else
			res.push_back(((t_series[i] - lmean[i]) * (1- (est_n/ float(lvar[i])))) + lmean[i]);
	}	
	return res;
}

int main(int argc, char *argv[]){
	std::vector<double> time_series = {112,31,42,54,432,245,23,34,567,8,468,623,232,743,121};
	int order = 2;
	std::vector<double> wf1d = wiener_filter_1d(time_series, order);
	std::cout.precision(12);
	for (auto it = wf1d.begin(); it != wf1d.end(); ++it)
		std::cout << *it << " ";
	std::cout << '\n';
}