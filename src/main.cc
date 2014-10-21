#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

typedef vector<vector<double> > d_mat;
#define INIT_DMAT(X,M,N,V) d_mat X(M,std::vector<double>(N,V))

vector<int> &char2Index(string const &s)
{
    static vector<int> n;
    int l = s.size();
    n.resize(l);
    for (int i =0; i < l; ++ i)
        n[i] = (int)s[i] - 65;
    return n;
}

double logsum(double const x, double const y) /* more strict log sum than naive */
{
    if (-INFINITY == y)
        return x;
    return x < y ? y + log1p(exp(x - y)) : x + log1p(exp(y - x));
}

    template <typename T>
void operator+= (vector<T> &v1, vector<T> const &v2)
{
    int l = v1.size() < v2.size()? v1.size(): v2.size();
    for (int i = 0; i < l; ++ i)
        v1[i] += v2[i];
}

    template <typename T1>
vector<T1> operator* (vector<T1> const &v, double const s)
{
    vector<T1> u(v.size());
    for (unsigned int i = 0; i < v.size(); ++ i)
        u[i] = v[i] * s;
    return u;
}

d_mat computeForward(vector<int> &s, d_mat model)
{
    int l = s.size(), mlen = model.size() - 1;
    INIT_DMAT(forward, l + 1, mlen + 2, -INFINITY);
    forward[0][0] = .0; /* initial state in the begining */
    for (int i = 1; i <= l; ++ i)
        for (int j = 0; j <= mlen + 1; ++ j)
            if (0 == j)
                forward[i][j] = forward[i - 1][j] + model[0][s[i - 1]];
            else if (1 <= j && j <= mlen)
                forward[i][j] = forward[i - 1][j - 1] + model[j][s[i - 1]];
            else /* mlen + 1 == j */
                forward[i][j] = logsum(forward[i - 1][j] + model[0][s[i - 1]], forward[i][j - 1]);
    return forward;
}

d_mat computeBackward(vector<int> &s, d_mat model)
{
    int l = s.size(), mlen = model.size() - 1;
    INIT_DMAT(backward, l + 1, mlen + 2, -INFINITY);
    backward[l] = vector<double>(mlen + 2, .0); /* any state in the end */
    for (int i = l - 1; i >= 0; -- i)
        for (int j = mlen + 1; j >= 0; -- j)
            if (mlen + 1 == j)
                backward[i][j] = model[0][s[i]] + backward[i + 1][j];
            else if (1 <= j && j <= mlen)
                backward[i][j] = model[j][s[i]] + backward[i + 1][j + 1];
            else /* 0 == j */
                backward[i][j] = logsum(model[0][s[i]] + backward[i + 1][j], backward[i][j + 1]);
    return backward;
}

void addExpCountMeme(vector<int> &s, d_mat &forward, d_mat &backward, d_mat &cv)
{
    int l = s.size(), mlen = cv.size() - 1, nchar = cv.front().size();
    vector<double> null0(nchar, .0);
    for (int i = 0; i < l; ++ i)
        ++ null0[s[i]]; /* base of null model */
    for (int i = 0; i <= l - mlen; ++ i) {
        double pr_m = exp(forward[i][1] + backward[i][1] - forward[l][mlen + 1]); /* prob that motif begins at i */
        vector<double> null = null0;
        for (int j = 0; j < mlen; ++ j)
            -- null[s[i + j]]; /* exclude inside of motif from null model */
        cv[0] += null * pr_m;
        for (int j = 0; j < mlen; ++ j)
            cv[j + 1][s[i + j]] += pr_m;
    }
}

d_mat computeExpectationMeme (vector<string> &seqs, d_mat &model)
{
    int mlen = model.size() - 2, nchar = model.front().size();
    INIT_DMAT(cv, mlen + 2, nchar, .01); /* initial value is as pseudocount */
    for (const auto &seq : seqs) { /* c++ 11 */
        vector<int> s = char2Index(seq);
        d_mat forward = computeForward(s, model);
        d_mat backward = computeBackward(s, model);
        addExpCountMeme(s, forward, backward, cv);
    }
    return cv;
}

d_mat solveEmMstep(d_mat &motif)
{
    int mlen = motif.size() - 1, nchar = motif.front().size();
    for (int i = 0; i < mlen; ++ i) {
        double z = .0;
        for (int c = 0; c < nchar; ++ c)
            z += motif[i][c];
        for (int c = 0; c < nchar; ++ c)
            motif[i][c] /= z;
    }
    return motif;
}

double dif(d_mat const &m1, d_mat const &m2)
{
    double z = .0;
    int m = m1.size(), n = m1.front().size();
    for (int i = 0; i < m; ++ i)
        for (int j = 0; j < n; ++ j)
            z += abs(m1[i][j] - m2[i][j]);
    return z;
}

d_mat trainEmMeme (vector<string> &seqs, d_mat model, double eps, int maxIt)
{
    while (maxIt --) { /* no limit (but overflow) with negative maxIt */
        d_mat cv = computeExpectationMeme(seqs, model);
        d_mat modelNew = solveEmMstep(cv);
        if (dif(model, modelNew) < eps)
            break;
        model = modelNew;
    }
    return model;
}

void printMostLikelyMotif(d_mat &motif)
{
    int mlen = motif.size() - 1, nchar = motif.front().size();
    for (int i = 1; i <= mlen; ++ i) {
        double max = 0;
        int maxIndex = 0;
        for (int j = 0; j < nchar; ++ j)
            if (max < motif[i][j])
                max = motif[i][maxIndex = j];
        cout << (char)(65 + maxIndex);
    }
    cout << endl;
}

int main(const int, const char *argv[])
{
    ifstream fasta(argv[1]);
    int mlen = atoi(argv[2]); /* motif length */
    int nchar = 26; /* how many kinds of chars */

    vector<string> seqs;
    for (string oneLine; getline(fasta, oneLine);)
        if ('>' == oneLine.at(0)) /* fasta header */
            seqs.push_back("");
        else
            seqs.back() += oneLine;

    INIT_DMAT(model, mlen + 1, nchar, 1.0 / nchar);
    d_mat motif = trainEmMeme(seqs, model, 0.01, 1000);

    printMostLikelyMotif(motif);

    return 0;
}
