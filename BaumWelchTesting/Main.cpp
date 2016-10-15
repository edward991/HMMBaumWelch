#include"../HMM/HMM.h"
#include <fstream>
#include <iostream>
#include <locale>
#include <codecvt>
using namespace std;

wchar_t *readObservationsFromFile(const char *file, int size)
{
	int i = 0;
	wchar_t *observations = new wchar_t[size];
	wifstream stream(file, ios::binary);

	stream.imbue(locale(stream.getloc(), new codecvt_utf8<wchar_t, 0x10ffff, consume_header>));

	for (wchar_t c; stream.get(c);)
	{
		if (c != ' ' && i < size)
		{
			observations[i] = c;
			i++;
		}
	}

	stream.close();

	return observations;
}

int numberOfWordsInFile(const char *file)
{
	int i = 0, j = 0;
	wifstream stream(file, ios::binary);

	stream.imbue(locale(stream.getloc(), new codecvt_utf8<wchar_t, 0x10ffff, consume_header>));

	for (wchar_t c; stream.get(c);)
	{
		if (c != ' ' && c != '\n' && c != '\r')
			j++;
		else
		{
			if (j > 0)
			{
				i++;
				j = 0;
			}
		}
	}

	stream.close();

	return i;
}

vector<vector<wchar_t>> readTrainingDataFromFile(const char *file)
{
	vector<vector<wchar_t>> trainingData(numberOfWordsInFile(file));
	int i = 0;
	vector<wchar_t> tempData;
	wifstream stream(file, ios::binary);

	stream.imbue(locale(stream.getloc(), new codecvt_utf8<wchar_t, 0x10ffff, consume_header>));

	for (wchar_t c; stream.get(c);)
	{
		if (c != ' ' && c != '\n' && c != '\r')
			tempData.push_back(c);
		else
		{
			if (tempData.size() > 0)
			{
				//trainingData.push_back(tempData);
				trainingData[i] = tempData;
				i++;
				tempData.clear();
			}
		}
	}

	stream.close();

	return trainingData;
}

template <class S, class O>
void printParametersToFile(const char *file, const HMM<S, O> &hmm)
{
	wofstream stream(file, ios::binary);

	stream.imbue(locale(stream.getloc(), new codecvt_utf8<wchar_t, 0x10ffff, consume_header>));

	stream << "Initial probabilities:\n";
	for (int i = 0; i < hmm.getNumberOfStates(); i++)
		stream << hmm.getInitialProbabilities()[i] << "  ";

	stream << "\n";

	stream << "\nTransition matrix:\n";
	for (int i = 0; i < hmm.getNumberOfStates(); i++)
	{
		for (int j = 0; j < hmm.getNumberOfStates(); j++)
			stream << hmm.getTransitionMatrix()[i][j] << "  ";
		stream << "\n";
	}

	stream << "\n";

	stream << "Emission matrix:\n";
	for (int i = 0; i < hmm.getNumberOfObservations(); i++)
	{
		stream << hmm.getObservations()[i] << " ";
		for (int j = 0; j < hmm.getNumberOfStates(); j++)
			stream << hmm.getEmissionMatrix()[j][i] << "  ";
		stream << '\n';
	}

	stream.close();
}

template <class S, class O>
void decodingTestData(const char *dataFile, const char *outputFile, HMM<S, O> &hmm)
{
	vector<vector<wchar_t>> data = readTrainingDataFromFile(dataFile);
	wchar_t *test;
	int *result;
	wofstream stream(outputFile, ios::binary);

	stream.imbue(locale(stream.getloc(), new codecvt_utf8<wchar_t, 0x10ffff, consume_header>));

	for (unsigned int i = 0; i < data.size(); i++)
	{
		test = new wchar_t[data[i].size()];
		result = new int[data[i].size()];

		for (unsigned int j = 0; j < data[i].size(); j++)
			test[j] = data[i][j];

		result = hmm.forwardBackwardAlgorithm(test, data[i].size());

		for (unsigned int j = 0; j < data[i].size(); j++)
			stream << test[j];
		stream << "\n";

		for (unsigned int j = 0; j < data[i].size(); j++)
			stream << result[j];
		stream << "\n";

		delete[] test;
		delete[] result;
	}

	stream.close();
}

int main()
{
	wchar_t *observations = readObservationsFromFile(".//Corpus/serbianCyrilicAlphabet.txt", 30);
	HMM<int, wchar_t> hmm = HMM<int, wchar_t>(2, observations, 30);

	vector<vector<wchar_t>> data = readTrainingDataFromFile(".//Corpus/trainingData.txt");

	hmm.trainParametersWithBaumWelchAlgorithm(data, 200, 0.000001);

	printParametersToFile<int, wchar_t>(".//Results/parameters.txt", hmm);
	decodingTestData<int, wchar_t>(".//Corpus/testData.txt", ".//Results/decodedData.txt", hmm);

	return 0;
}