//
// Created by rbd on 7.6.2020.
//

#ifndef EPSILONNET_RUN1MEANCORE_H
#define EPSILONNET_RUN1MEANCORE_H

#include <vector>
#include <string>

using namespace std;

void toCSV(vector<vector<double> > P, int n, int d, string file);

vector<vector<double> > & parseCSV(string file);

//inline vector<vector<double> > & runCoreset(
 vector<vector<double> > & runCoreset(
        vector<vector<double> > & P, int n, int d, double eps,
        double alpha = 1, double delta = 0.1,
        bool isPrivate = true, int security = 1024);


#endif //EPSILONNET_RUN1MEANCORE_H
