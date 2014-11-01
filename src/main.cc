#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
typedef std::vector<std::vector<double> > d_mat;

    template <typename T>
void operator+= (std::vector<T> &v1, std::vector<T> const &v2)
{
    int l = v1.size() < v2.size()? v1.size(): v2.size();
    for (int i = 0; i < l; ++ i)
        v1[i] += v2[i];
}
    template <typename T>
std::vector<T> operator* (std::vector<T> const &v, double const s)
{
    std::vector<T> u(v.size());
    for (size_t i = 0; i < (int)v.size(); ++ i)
        u[i] = v[i] * s;
    return u;
}

struct Meme
{
    std::vector<std::string> seqs; /* sequences */
    std::vector<std::vector<int>> nseqs; /* integer form sequence */
    std::string const alphs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    int mlen, nalph = alphs.size();/* length of motif, number of alphabets */
    d_mat model, forward, backward;
    Meme(std::vector<std::string> const &s, int const ml): seqs(s), mlen(ml)
    {
        nseqs.resize(s.size());
        for (size_t i = 0; i < s.size(); ++ i) {
            nseqs[i].resize(seqs[i].size());
            std::transform(s[i].begin(), s[i].end(), nseqs[i].begin(), [&] (char c) {return alphs.find(c);});
        }
        model.resize(ml + 1);
        std::fill_n(model.begin(), ml + 1, std::vector<double>(nalph, 1.0 / nalph));
    }
    double logsum(double const x, double const y) /* more strict log(sum(exp..)) than naive */
    {
        if (-INFINITY == y)
            return x;
        return x < y ? y + std::log1p(exp(x - y)) : x + std::log1p(exp(y - x));
    }
    void computeForwardBackward(std::vector<int> &s)
    {
        forward.resize(s.size() + 1);
        std::fill_n(forward.begin(), s.size() + 1, std::vector<double>(mlen + 2, -INFINITY));
        forward[0][0] = .0; /* initial state in the begining */
        for (size_t i = 1; i <= s.size(); ++ i)
            for (int j = 0; j <= mlen + 1; ++ j)
                if (0 == j)
                    forward[i][j] = forward[i - 1][j] + log(model[0][s[i - 1]]);
                else if (1 <= j && j <= mlen)
                    forward[i][j] = forward[i - 1][j - 1] + log(model[j][s[i - 1]]);
                else /* mlen + 1 == j */
                    forward[i][j] = logsum(forward[i - 1][j], forward[i - 1][j - 1]) + log(model[0][s[i - 1]]);
        backward.resize(s.size() + 1);
        std::fill_n(backward.begin(), s.size() + 1, std::vector<double>(mlen + 2, -INFINITY));
        backward[s.size()] = std::vector<double>(mlen + 2, .0); /* any state in the end */
        for (int i = s.size() - 1; i >= 0; -- i)
            for (int j = mlen + 1; j >= 0; -- j)
                if (mlen <= j)
                    backward[i][j] = backward[i + 1][mlen + 1] + log(model[0][s[i]]);
                else if (1 <= j && j < mlen)
                    backward[i][j] = backward[i + 1][j + 1] + log(model[j + 1][s[i]]);
                else /* 0 == j */
                    backward[i][j] = logsum(backward[i + 1][j] + log(model[0][s[i]]), backward[i + 1][j + 1] + log(model[j + 1][s[i]]));
    }
    void addExpCountMeme(std::vector<int> const &s, d_mat &cv)
    {
        std::vector<double> null0(nalph, .0);
        for (size_t i = 0; i < s.size(); ++ i)
            null0[s[i]] += 1.0; /* base of null model */
        for (size_t i = 0; i <= s.size() - mlen; ++ i) {
            double pr_m = exp(forward[i + 1][1] + backward[i + 1][1] - forward[s.size()][mlen + 1]); /* prob that motif begins at i */
            std::vector<double> null = null0;
            for (int j = 0; j < mlen; ++ j)
                null[s[i + j]] -= 1.0; /* exclude inside of motif from null model */
            cv[0] += null * pr_m;
            for (int j = 0; j < mlen; ++ j)
                cv[j + 1][s[i + j]] += pr_m;
        }
    }
    d_mat computeExpectation(double const pseudocount)
    {
        d_mat cv(mlen + 1, std::vector<double>(nalph, pseudocount)); /* initial value is as pseudocount */
        for (size_t i = 0; i < seqs.size(); ++ i) {
            computeForwardBackward(nseqs[i]);
            addExpCountMeme(nseqs[i], cv);
        }
        return cv;
    }
    d_mat solveEmMstep(d_mat &motif)
    {
        for (int i = 0; i <= mlen; ++ i) {
            double z = .0;
            for (int c = 0; c < nalph; ++ c)
                z += motif[i][c];
            for (int c = 0; c < nalph; ++ c)
                motif[i][c] /= z;
        }
        return motif;
    }
    double dif(d_mat const &m1, d_mat const &m2)
    {
        double z = .0;
        for (size_t i = 0; i < m1.size(); ++ i)
            for (size_t j = 0; j < m1.front().size(); ++ j)
                z += std::abs(m1[i][j] - m2[i][j]);
        return z;
    }
    void trainEm(double const eps, int maxIt, double const pseudocount)
    {
        while (maxIt --) { /* no limit (but overflow) with negative maxIt */
            d_mat cv = computeExpectation(pseudocount);
            d_mat modelNew = solveEmMstep(cv);
            if (dif(model, modelNew) < eps)
                break;
            model = modelNew;
        }
    }
    std::string mostLikelyMotif(void)
    {
        std::string motif_seq;
        for (int i = 1; i <= mlen; ++ i) {
            double max = .0;
            int maxIndex = 0;
            for (int j = 0; j < nalph; ++ j)
                if (max < model[i][j])
                    max = model[i][maxIndex = j];
            motif_seq += alphs[maxIndex];
        }
        return motif_seq;
    }
};

std::vector<std::string> parseFasta(std::ifstream &fasta)
{
    std::vector<std::string> seqs;
    for (std::string oneLine; getline(fasta, oneLine);)
        if ('>' == oneLine.at(0)) /* fasta header */
            seqs.push_back("");
        else
            seqs.back() += oneLine;
    return seqs;
}
int main(const int, const char *argv[])
{
    std::ifstream fasta(argv[1]);
    std::vector<std::string> seqs = parseFasta(fasta);
    Meme meme(seqs, atoi(argv[2]));
    meme.trainEm(.01, 1000, .01);
    std::cout << meme.mostLikelyMotif() << std::endl;
    return 0;
}
