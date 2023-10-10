
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct constant_group {
    int random_seed_genotype = 0;
    int random_seed_antibiotics = 0;
    double alpha = 0.0;
    int CF_index = 0;
    int index_cell = 0;
    int Import_Max_Evolution_Time = 5000;
    int Rank_lowdim = 0;
    int N_gene = 100;
    int N_edge = 100;
    int N_input = 8;
    int N_target = 8;
    double beta_s = 1.0;
    double beta_t = 1.0;
    double theta = 0.0;
    double s_I = 1.0 / N_input;
    vector<vector<double>> Input;
    vector<vector<double>> Target;
    vector<double> Antibiotics;
    int N_cell = 100;
    int N_selection = 25;
    int Max_Evolution_Time = 2000;
    double N_mutation = 1.0;

    // Parameters of rungekutta
    double initial_dt = 0.01;
    double max_dt = 0.1;
    double min_dt = 1.0 / pow(10, 8);
    double min_error = 1.0 / pow(10, 4);
    double tol_error = 1.0 / pow(10, 4);
    double a21 = 1 / 4.0;
    double a31 = 3.0 / 32.0;
    double a32 = 9.0 / 32.0;
    double a41 = 1932.0 / 2197.0;
    double a42 = -7200.0 / 2197.0;
    double a43 = 7296.0 / 2197.0;
    double a51 = 439.0 / 216.0;
    double a52 = -8.0;
    double a53 = 3680.0 / 513.0;
    double a54 = -845.0 / 4104.0;
    double a61 = -8.0 / 27.0;
    double a62 = 2.0;
    double a63 = -3544.0 / 2565.0;
    double a64 = 1859.0 / 4104.0;
    double a65 = -11.0 / 40.0;
};

double Calc_fitness(const constant_group& cg, const vector<vector<int>>& Genome_cell, vector<double>& x_cell);
int adaptive_rungekutta(const constant_group& cg, vector<double>& x, const vector<vector<int>>& Genome_cell, double& dt, double& t);
void calc_k(const constant_group& cg, vector<double>& x, const vector<vector<int>>& Genome_cell, double& dt, vector<vector<double>>& k, int rank_k);
void make_x_dummy(const constant_group cg, vector<double>& x_dummy, vector<double>& x, vector<vector<double>>& k, int rank_k);
int pow_int(int a, int b);
void Make_Antibiotics(constant_group& cg);
int import_Genome(vector<vector<int>>& Genome_cell, string filename);

double Calc_fitness(const constant_group& cg, const vector<vector<int>>& Genome_cell, vector<double>& x_cell)
{
    random_device rnd;
    mt19937_64 mt(rnd());
    uniform_real_distribution<> rnd_01(0.0, 1.0);
    double t = 0.0, dt = cg.initial_dt;
    double tmp_target = 0.0;
    vector<double> fitness_vector(cg.N_target, 0.0);
    double x_tmp = 0.0;
    int index_t = 0;
    for (int i = 0; i < cg.N_gene; ++i) {
        x_cell[i] = rnd_01(mt);
    }
    while (adaptive_rungekutta(cg, x_cell, Genome_cell, dt, t) == 0) {
    }
    for (int i = 0; i < cg.N_target; ++i) {
        x_tmp = 0.0;
        for (int j = 0; j < cg.N_gene; ++j) {
            x_tmp += Genome_cell[cg.N_gene + i][j] * x_cell[j] / sqrt((double)cg.N_edge);
        }
        x_tmp = 1.0 / (1.0 + exp(-cg.beta_t * x_tmp));
        fitness_vector[i] = (x_tmp - cg.Target[cg.CF_index][i]) * (x_tmp - cg.Target[cg.CF_index][i]);
    }
    return accumulate(fitness_vector.begin(), fitness_vector.end(), 0.0);
}
int adaptive_rungekutta(const constant_group& cg, vector<double>& x, const vector<vector<int>>& Genome_cell, double& dt, double& t)
{
    vector<vector<double>> k(6, vector<double>(cg.N_gene, 0.0));
    for (int i = 0; i < 6; ++i) {
        calc_k(cg, x, Genome_cell, dt, k, i);
    }
    double k2_sum = 0.0;
    for (int i = 0; i < cg.N_gene; ++i) {
        k2_sum += (25.0 / 216.0 * k[0][i] + 1408.0 / 2565.0 * k[2][i] + 2197.0 / 4104.0 * k[3][i] - 0.20 * k[4][i]) * (25.0 / 216.0 * k[0][i] + 1408.0 / 2565.0 * k[2][i] + 2197.0 / 4104.0 * k[3][i] - 0.20 * k[4][i]);
    }
    k2_sum = sqrt(k2_sum) / dt;
    if (k2_sum < cg.min_error) {
        return 1;
    }
    double R = 0.0;
    for (int i = 0; i < cg.N_gene; ++i) {
        R += (k[0][i] / 360.0 - k[2][i] * 128.0 / 4275.0 - k[3][i] * 2197.0 / 75240.0 + k[4][i] / 50.0 + k[5][i] * 2.0 / 55.0) * (k[0][i] / 360.0 - k[2][i] * 128.0 / 4275.0 - k[3][i] * 2197.0 / 75240.0 + k[4][i] / 50.0 + k[5][i] * 2.0 / 55.0);
    }
    R = sqrt(R) / dt;
    if (R <= cg.tol_error) {
        t += dt;
        double sum_x = 0.0;
        for (int i = 0; i < cg.N_gene; ++i) {
            x[i] += 25.0 / 216.0 * k[0][i];
            x[i] += 1408.0 / 2565.0 * k[2][i];
            x[i] += 2197.0 / 4104.0 * k[3][i];
            x[i] -= 0.20 * k[4][i];
        }
    }
    double delta_x = pow(cg.tol_error / 2.0 / R, 0.25);
    if (delta_x <= 0.1) {
        dt *= 0.1;
    } else if (delta_x >= 4.0) {
        dt *= 4.0;
    } else {
        dt *= delta_x;
    }
    if (dt > cg.max_dt) {
        dt = cg.max_dt;
    } else if (dt < cg.min_dt) {
        dt = cg.min_dt;
    }
    return 0;
}
void calc_k(const constant_group& cg, vector<double>& x, const vector<vector<int>>& Genome_cell, double& dt, vector<vector<double>>& k, int rank_k)
{
    vector<double> x_dummy(cg.N_gene, 0.0);
    make_x_dummy(cg, x_dummy, x, k, rank_k);
    vector<double> x_tilde(cg.N_gene, 0.0);
    for (int i = 0; i < cg.N_gene; ++i) {
        for (int j = 0; j < cg.N_gene; ++j) {
            x_tilde[i] += Genome_cell[i][j] * x_dummy[j] / sqrt((double)cg.N_edge);
        }
        x_tilde[i] += cg.Antibiotics[i];
        k[rank_k][i] = (1.0 / (1.0 + exp(-cg.beta_s * x_tilde[i])) - x_dummy[i]) * dt;
    }
    return;
}
void make_x_dummy(const constant_group cg, vector<double>& x_dummy, vector<double>& x, vector<vector<double>>& k, int rank_k)
{
    if (rank_k == 0) {
        x_dummy = x;
    } else if (rank_k == 1) {
        for (int i = 0; i < cg.N_gene; ++i) {
            x_dummy[i] = x[i] + cg.a21 * k[0][i];
        }
    } else if (rank_k == 2) {
        for (int i = 0; i < cg.N_gene; ++i) {
            x_dummy[i] = x[i] + cg.a31 * k[0][i] + cg.a32 * k[1][i];
        }
    } else if (rank_k == 3) {
        for (int i = 0; i < cg.N_gene; ++i) {
            x_dummy[i] = x[i] + cg.a41 * k[0][i] + cg.a42 * k[1][i] + cg.a43 * k[2][i];
        }
    } else if (rank_k == 4) {
        for (int i = 0; i < cg.N_gene; ++i) {
            x_dummy[i] = x[i] + cg.a51 * k[0][i] + cg.a52 * k[1][i] + cg.a53 * k[2][i] + cg.a54 * k[3][i];
        }
    } else if (rank_k == 5) {
        for (int i = 0; i < cg.N_gene; ++i) {
            x_dummy[i] = x[i] + cg.a61 * k[0][i] + cg.a62 * k[1][i] + cg.a63 * k[2][i] + cg.a64 * k[3][i] + cg.a65 * k[4][i];
        }
    }
    return;
}
int pow_int(int a, int b)
{
    int x = 1;
    for (int i = 0; i < b; ++i) {
        x *= a;
    }
    return x;
}
void Make_Antibiotics(constant_group& cg)
{
    mt19937_64 mt(cg.random_seed_antibiotics);
    normal_distribution<> Gauss(0.0, 1.0);
    vector<double> tmp_antibiotics(cg.N_gene, 0.0);
    for (int i = 0; i < cg.N_gene; ++i) {
        tmp_antibiotics[i] = Gauss(mt);
    }
    cg.Antibiotics = tmp_antibiotics;
    return;
}
int import_Genome(vector<vector<int>>& Genome_cell, string filename)
{
    vector<vector<int>> edge_temp;
    ifstream ifs(filename);
    if (!ifs) {
        cout << "入力エラー" << endl;
        return 0;
    }
    string str;
    while (getline(ifs, str)) {
        string token;
        istringstream stream(str);
        vector<int> tmp;
        while (getline(stream, token, ',')) {
            tmp.push_back(stoi(token));
        }
        edge_temp.push_back(tmp);
    }
    Genome_cell = edge_temp;
    return 1;
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        cout << "引数の数がおかしい" << endl;
        return 0;
    }
    constant_group cg;
    cg.Rank_lowdim = atoi(argv[1]);
    cg.alpha = atof(argv[2]);
    cg.random_seed_genotype = atoi(argv[3]);
    cg.index_cell = atoi(argv[4]);

    vector<vector<vector<int>>> Genome(cg.N_cell, vector<vector<int>>(cg.N_gene + cg.N_target, vector<int>(cg.N_gene + cg.N_input, 0)));
    vector<vector<vector<int>>> Genome_selection(cg.N_selection, vector<vector<int>>(cg.N_gene + cg.N_target, vector<int>(cg.N_gene + cg.N_input, 0)));
    vector<vector<double>> x(cg.N_cell, vector<double>(cg.N_gene, 0.0));
    vector<pair<double, int>> fitness(cg.N_cell);

    if (import_Genome(Genome_selection[cg.index_cell], "./Np_" + to_string(cg.Rank_lowdim) + "_" + to_string(cg.alpha) + "_" + to_string(cg.random_seed_genotype) + "/Genome_" + to_string(cg.index_cell) + "_.csv") == 0) {
        return 0;
    }
    ofstream outputfile_fitness("./Np_" + to_string(cg.Rank_lowdim) + "_" + to_string(cg.alpha) + "_" + to_string(cg.random_seed_genotype) + "/VE_fitness_" + to_string(cg.index_cell) + "_.csv");
    ofstream outputfile_x("./Np_" + to_string(cg.Rank_lowdim) + "_" + to_string(cg.alpha) + "_" + to_string(cg.random_seed_genotype) + "/VE_x_" + to_string(cg.index_cell) + "_.csv");

    cg.Input.emplace_back(vector<double>(cg.N_input, 0.0));
    cg.Target.emplace_back(vector<double>(cg.N_target, 0.5));
    for (cg.random_seed_antibiotics = 0; cg.random_seed_antibiotics < 10000; cg.random_seed_antibiotics += 1) {
        Make_Antibiotics(cg);
        outputfile_fitness << Calc_fitness(cg, Genome_selection[cg.index_cell], x[cg.index_cell]) << endl;
        outputfile_x << x[cg.index_cell][0];
        for (int j = 1; j < cg.N_gene; ++j) {
            outputfile_x << "," << x[cg.index_cell][j];
        }
        outputfile_x << endl;
    }
    return 0;
}