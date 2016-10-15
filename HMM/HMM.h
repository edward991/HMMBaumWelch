#pragma once
#include<iostream>
#include<vector>
#include<cstdlib>
#include<cmath>
using namespace std;

template <class S, class O>
class HMM
{
private:
	static const int MAX_OBS_SIZE = 100;

	int numberOfStates;
	int numberOfObservations;

	S *states;
	O *observations;
	double **transitionMatrix;
	double **emissionMatrix;
	double *initialProbabilities;

	double **alpha;
	double **beta;
	double **gamma;
	double ***epsilon;

public:

	HMM(int numberOfStates, O *observations, int numberOfObservations)
	{
		this->numberOfStates = numberOfStates;
		this->numberOfObservations = numberOfObservations;

		states = new S[numberOfStates];
		for (int i = 0; i < numberOfStates; i++)
			states[i] = (S)i;

		this->observations = observations;

		transitionMatrix = new double*[numberOfStates];
		for (int i = 0; i < numberOfStates; i++)
			transitionMatrix[i] = new double[numberOfStates];

		emissionMatrix = new double*[numberOfStates];
		for (int i = 0; i < numberOfStates; i++)
			emissionMatrix[i] = new double[numberOfObservations];

		initialProbabilities = new double[numberOfStates];

		this->initialize();
	}

	HMM(S *states, int numberOfStates, O *observations, int numberOfObservations, double **transitionMatrix, double **emissionMatrix, double *initialProbabilities)
	{
		this->numberOfStates = numberOfStates;
		this->numberOfObservations = numberOfObservations;

		this->states = states;
		this->observations = observations;
		this->transitionMatrix = transitionMatrix;
		this->emissionMatrix = emissionMatrix;
		this->initialProbabilities = initialProbabilities;

		this->initialize();
	}

	~HMM()
	{
		delete[] this->states;
		delete[] this->observations;
		delete[] this->initialProbabilities;

		for (int i = 0; i < numberOfStates; i++)
		{
			delete[] this->transitionMatrix[i];
			delete[] this->emissionMatrix[i];

			for (int j = 0; j < numberOfStates; j++)
				delete[] this->epsilon[i][j];

			delete[] this->epsilon[i];
		}

		for (int i = 0; i < MAX_OBS_SIZE; i++)
		{
			delete[] this->alpha[i];
			delete[] this->beta[i];
			delete[] this->gamma[i];
		}

		delete[] this->transitionMatrix;
		delete[] this->emissionMatrix;
		delete[] this->alpha;
		delete[] this->beta;
		delete[] this->gamma;
		delete[] this->epsilon;
	}

	int getNumberOfStates() const
	{
		return this->numberOfStates;
	}

	int getNumberOfObservations() const
	{
		return this->numberOfObservations;
	}

	S* getStates() const
	{
		return this->states;
	}

	O* getObservations() const
	{
		return this->observations;
	}

	double** getTransitionMatrix() const
	{
		return this->transitionMatrix;
	}

	double** getEmissionMatrix() const
	{
		return this->emissionMatrix;
	}

	double* getInitialProbabilities() const
	{
		return this->initialProbabilities;
	}

private:

	void initialize()
	{
		alpha = new double*[MAX_OBS_SIZE];
		for (int i = 0; i < MAX_OBS_SIZE; i++)
			alpha[i] = new double[numberOfStates];

		beta = new double*[MAX_OBS_SIZE];
		for (int i = 0; i < MAX_OBS_SIZE; i++)
			beta[i] = new double[numberOfStates];

		gamma = new double*[MAX_OBS_SIZE];
		for (int i = 0; i < MAX_OBS_SIZE; i++)
			gamma[i] = new double[numberOfStates];

		epsilon = new double**[numberOfStates];
		for (int i = 0; i < numberOfStates; i++)
		{
			epsilon[i] = new double*[numberOfStates];
			for (int j = 0; j < numberOfStates; j++)
				epsilon[i][j] = new double[MAX_OBS_SIZE];
		}
	}

	int positionInObservations(O o) const
	{
		for (int i = 0; i < numberOfObservations; i++)
		{
			if (observations[i] == o)
				return i;
		}

		return -1;
	}

	int positionInStates(S s) const
	{
		for (int i = 0; i < numberOfStates; i++)
		{
			if (states[i] == s)
				return i;
		}

		return -1;
	}

	double* generateRandomDistribution(int size, int range) const
	{
		double *randomArray = new double[size];
		int randomNumber, sum = 0;

		for (int i = 0; i < size; i++)
		{
			randomNumber = rand() % range + 1;
			sum += randomNumber;
			randomArray[i] = (double)randomNumber;
		}

		for (int i = 0; i < size; i++)
			randomArray[i] /= (double)sum;

		return randomArray;
	}

	void generateRandomParameters(int range)
	{
		delete[] initialProbabilities;
		initialProbabilities = generateRandomDistribution(numberOfStates, range);

		for (int i = 0; i < numberOfStates; i++)
		{
			delete[] transitionMatrix[i];
			transitionMatrix[i] = generateRandomDistribution(numberOfStates, range);
		}

		for (int i = 0; i < numberOfStates; i++)
		{
			delete[] emissionMatrix[i];
			emissionMatrix[i] = generateRandomDistribution(numberOfObservations, range);
		}
	}

	void evaluateAlpha(const O *observationSequence, int t)
	{
		double value;

		int *positions = new int[t];
		for (int i = 0; i < t; i++)
			positions[i] = positionInObservations(observationSequence[i]);

		for (int k = 0; k < t; k++)
		{
			for (int i = 0; i < numberOfStates; i++)
			{
				value = 0.0;
				if (k == 0)
					value = initialProbabilities[i] * emissionMatrix[i][positions[0]];
				else
				{
					for (int j = 0; j < numberOfStates; j++)
						value += alpha[k - 1][j] * transitionMatrix[j][i] * emissionMatrix[i][positions[k]];
				}
				alpha[k][i] = value;
			}
		}

		delete[] positions;
	}

	void evaluateBeta(const O *observationSequence, int t)
	{
		double value;

		int *positions = new int[t];
		for (int i = 0; i < t; i++)
			positions[i] = positionInObservations(observationSequence[i]);

		for (int k = t - 1; k >= 0; k--)
		{
			for (int i = 0; i < numberOfStates; i++)
			{
				value = 0.0;
				if (k == t - 1)
					value = 1.0;
				else
				{
					for (int j = 0; j < numberOfStates; j++)
						value += beta[k + 1][j] * transitionMatrix[i][j] * emissionMatrix[j][positions[k + 1]];
				}
				beta[k][i] = value;
			}
		}

		delete[] positions;
	}

	void evaluateGamma(const O *observationSequence, int t)
	{
		double value, sum;

		sum = 0.0;
		for (int i = 0; i < numberOfStates; i++)
			sum += alpha[t - 1][i];

		for (int k = 0; k < t; k++)
		{
			for (int i = 0; i < numberOfStates; i++)
			{
				value = alpha[k][i] * beta[k][i];
				gamma[k][i] = value / sum;
			}
		}
	}

	void evaluateEpsilon(const O *observationSequence, int t)
	{
		double value, sum;
		sum = 0.0;

		int *positions = new int[t];
		for (int i = 0; i < t; i++)
			positions[i] = positionInObservations(observationSequence[i]);

		for (int i = 0; i < numberOfStates; i++)
			sum += alpha[t - 1][i];

		for (int i = 0; i < numberOfStates; i++)
		{
			for (int j = 0; j < numberOfStates; j++)
			{
				for (int k = 0; k < t - 1; k++)
				{
					value = alpha[k][i] * transitionMatrix[i][j] * emissionMatrix[j][positions[k + 1]] * beta[k + 1][j];
					epsilon[i][j][k] = value / sum;
				}
			}
		}

		delete[] positions;
	}

public:

	double forwardAlgorithm(const O *observationSequence, int t)
	{
		double p = 0.0;

		if (t > MAX_OBS_SIZE)
			t = MAX_OBS_SIZE;

		evaluateAlpha(observationSequence, t);

		for (int i = 0; i < numberOfStates; i++)
			p += alpha[t - 1][i];

		return p;
	}

	double backwardAlgorithm(const O *observationSequence, int t)
	{
		double p = 0.0;

		if (t > MAX_OBS_SIZE)
			t = MAX_OBS_SIZE;

		evaluateBeta(observationSequence, t);

		for (int i = 0; i < numberOfStates; i++)
			p += beta[0][i] * initialProbabilities[i] * emissionMatrix[i][positionInObservations(observationSequence[0])];

		return p;
	}

	S* forwardBackwardAlgorithm(const O *observationSequence, int t)
	{
		double max;
		S maxState;
		S *optimalStates = new S[t];

		if (t > MAX_OBS_SIZE)
			t = MAX_OBS_SIZE;

		evaluateAlpha(observationSequence, t);
		evaluateBeta(observationSequence, t);
		evaluateGamma(observationSequence, t);

		for (int k = 0; k < t; k++)
		{
			max = 0.0;

			for (int i = 0; i < numberOfStates; i++)
			{
				if (gamma[k][i] >= max)
				{
					max = gamma[k][i];
					maxState = states[i];
				}
			}

			optimalStates[k] = maxState;
		}
		return optimalStates;
	}

	S* viterbiAlgorithm(const O *observationSequence, int t)
	{
		double max, current;
		S maxState;

		double **delta = new double*[t];
		for (int i = 0; i < t; i++)
			delta[i] = new double[numberOfStates];

		S **omega = new S*[t];
		for (int i = 0; i < t; i++)
			omega[i] = new S[numberOfStates];

		S *optimalStates = new S[t];

		if (t > MAX_OBS_SIZE)
			t = MAX_OBS_SIZE;

		for (int i = 0; i < numberOfStates; i++)
		{
			delta[0][i] = initialProbabilities[i] * emissionMatrix[i][positionInObservations(observationSequence[0])];;
			//omega[0][i] = S();
		}

		for (int k = 1; k < t; k++)
		{
			for (int j = 0; j < numberOfStates; j++)
			{
				max = 0.0;

				for (int i = 0; i < numberOfStates; i++)
				{
					current = delta[k - 1][i] * transitionMatrix[i][j] * emissionMatrix[j][positionInObservations(observationSequence[k])];

					if (current >= max)
					{
						max = current;
						maxState = states[i];
					}
				}

				delta[k][j] = max;
				omega[k][j] = maxState;
			}
		}

		max = 0.0;

		for (int i = 0; i < numberOfStates; i++)
		{
			if (delta[t - 1][i] >= max)
			{
				max = delta[t - 1][i];
				optimalStates[t - 1] = states[i];
			}
		}

		for (int k = t - 1; k >= 1; k--)
			optimalStates[k - 1] = omega[k][positionInStates(optimalStates[k])];

		for (int i = 0; i < t; i++)
		{
			delete[] delta[i];
			delete[] omega[i];
		}

		delete[] delta;
		delete[] omega;

		return optimalStates;
	}

	void trainParametersByCountingFrequencies(const vector<vector<pair<O, S>>> &trainingData)
	{

		int count, count1, count2, numOfObservations = trainingData.size();

		for (int i = 0; i < numberOfStates; i++)
		{
			count = 0;

			for (int k = 0; k < numOfObservations; k++)
			{
				if (trainingData[k][0].second == states[i])
					count++;
			}

			initialProbabilities[i] = ((double)count / (double)numOfObservations);
		}

		for (int i = 0; i < numberOfStates; i++)
		{
			for (int j = 0; j < numberOfStates; j++)
			{
				count1 = 0;
				count2 = 0;

				for (int k = 0; k < numOfObservations; k++)
				{
					for (int l = 0; l < trainingData[k].size() - 1; l++)
					{
						if (trainingData[k][l].second == states[i])
						{
							count2++;
							if (trainingData[k][l + 1].second == states[j])
								count1++;
						}
					}
				}

				transitionMatrix[i][j] = (count2 == 0) ? 0 : ((double)count1 / (double)count2);
			}
		}

		for (int i = 0; i < numberOfStates; i++)
		{
			for (int j = 0; j < numberOfObservations; j++)
			{
				count1 = 0;
				count2 = 0;

				for (int k = 0; k < numOfObservations; k++)
				{
					for (int l = 0; l < trainingData[k].size(); l++)
					{
						if (trainingData[k][l].second == states[i])
						{
							count2++;
							if (trainingData[k][l].first == observations[j])
								count1++;
						}
					}
				}
				emissionMatrix[i][j] = (count2 == 0) ? 0 : ((double)count1 / (double)count2);
			}
		}
	}

	void trainParametersWithBaumWelchAlgorithm(const vector<vector<O>> &trainingData, int numOfIterations, double error)
	{
		double sum1, sum2, newLikelihood = 0.0, oldlikelihood = -10000;
		int iterations = 0, observationSequenceSize, numOfObservationSequences = trainingData.size();
		O *currentObservationSequnce;

		double *tempInitialProbability = new double[numberOfStates];

		double **tempTransitionMatrixNumerator = new double*[numberOfStates];
		double **tempTransitionMatrixDenominator = new double*[numberOfStates];
		for (int i = 0; i < numberOfStates; i++)
		{
			tempTransitionMatrixNumerator[i] = new double[numberOfStates];
			tempTransitionMatrixDenominator[i] = new double[numberOfStates];
		}

		double **tempEmissionMatrixNumerator = new double*[numberOfStates];
		double **tempEmissionMatrixDenominator = new double*[numberOfStates];
		for (int i = 0; i < numberOfStates; i++)
		{
			tempEmissionMatrixNumerator[i] = new double[numberOfObservations];
			tempEmissionMatrixDenominator[i] = new double[numberOfObservations];
		}

		generateRandomParameters(3);

		cout << "Iteration  |Log likelihood\n\n";
		while (true)
		{
			if (iterations == numOfIterations)
				break;
			iterations++;

			if (abs(newLikelihood - oldlikelihood) <= error)
				break;

			for (int i = 0; i < numberOfStates; i++)
				tempInitialProbability[i] = 0.0;

			for (int i = 0; i < numberOfStates; i++)
				for (int j = 0; j < numberOfStates; j++)
					tempTransitionMatrixNumerator[i][j] = tempTransitionMatrixDenominator[i][j] = 0.0;

			for (int i = 0; i < numberOfStates; i++)
				for (int j = 0; j < numberOfObservations; j++)
					tempEmissionMatrixNumerator[i][j] = tempEmissionMatrixDenominator[i][j] = 0.0;

			oldlikelihood = newLikelihood;
			newLikelihood = 0.0;

			for (int k = 0; k < numOfObservationSequences; k++)
			{
				observationSequenceSize = trainingData[k].size();
				if (observationSequenceSize > MAX_OBS_SIZE)
					observationSequenceSize = MAX_OBS_SIZE;

				currentObservationSequnce = new O[observationSequenceSize];
				for (int i = 0; i < observationSequenceSize; i++)
					currentObservationSequnce[i] = trainingData[k][i];

				evaluateAlpha(currentObservationSequnce, observationSequenceSize);
				evaluateBeta(currentObservationSequnce, observationSequenceSize);
				evaluateGamma(currentObservationSequnce, observationSequenceSize);
				evaluateEpsilon(currentObservationSequnce, observationSequenceSize);

				newLikelihood += log(forwardAlgorithm(currentObservationSequnce, observationSequenceSize));

				for (int i = 0; i < numberOfStates; i++)
				{
					tempInitialProbability[i] += gamma[0][i];

					for (int j = 0; j < numberOfStates; j++)
					{
						sum1 = 0.0;
						sum2 = 0.0;
						for (int l = 0; l < observationSequenceSize - 1; l++)
						{
							sum1 += epsilon[i][j][l];
							sum2 += gamma[l][i];
						}

						tempTransitionMatrixNumerator[i][j] += sum1;
						tempTransitionMatrixDenominator[i][j] += sum2;
					}

					for (int j = 0; j < numberOfObservations; j++)
					{
						sum1 = 0.0;
						sum2 = 0.0;
						for (int l = 0; l < observationSequenceSize; l++)
						{
							if (currentObservationSequnce[l] == observations[j])
								sum1 += gamma[l][i];
							sum2 += gamma[l][i];
						}
						tempEmissionMatrixNumerator[i][j] += sum1;
						tempEmissionMatrixDenominator[i][j] += sum2;
					}
				}

				delete[] currentObservationSequnce;
			}

			newLikelihood /= trainingData.size();

			cout << iterations << "          " << newLikelihood << "\n";

			for (int i = 0; i < numberOfStates; i++)
				initialProbabilities[i] = tempInitialProbability[i] / trainingData.size();

			for (int i = 0; i < numberOfStates; i++)
				for (int j = 0; j < numberOfStates; j++)
					transitionMatrix[i][j] = tempTransitionMatrixNumerator[i][j] / tempTransitionMatrixDenominator[i][j];

			for (int i = 0; i < numberOfStates; i++)
				for (int j = 0; j < numberOfObservations; j++)
					emissionMatrix[i][j] = tempEmissionMatrixNumerator[i][j] / tempEmissionMatrixDenominator[i][j];
		}

		delete[] tempInitialProbability;

		for (int i = 0; i < numberOfStates; i++)
		{
			delete[] tempTransitionMatrixNumerator[i];
			delete[] tempTransitionMatrixDenominator[i];
			delete[] tempEmissionMatrixNumerator[i];
			delete[] tempEmissionMatrixDenominator[i];
		}

		delete[] tempTransitionMatrixNumerator;
		delete[] tempTransitionMatrixDenominator;
		delete[] tempEmissionMatrixNumerator;
		delete[] tempEmissionMatrixDenominator;
	}
};

