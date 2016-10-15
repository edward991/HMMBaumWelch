#include"../HMM/HMM.h"

ostream &operator<<(ostream &stream, const string &s)
{
	int len = s.length();
	for (int i = 0; i < len; i++)
		stream << s[i];
	stream << " ";
	return stream;
}

int main()
{
	string *states=new string[2];
	states[0] = "Healthy";
	states[1] = "Fever";

	string *observations=new string[3];
	observations[0]="normal";
	observations[1]="cold";
	observations[2]="dizzy";

	double **transitionMatrix=new double*[2];
	for (int i = 0; i < 2; i++)
		transitionMatrix[i] = new double[2];
	transitionMatrix[0][0] = 0.7;
	transitionMatrix[0][1] = 0.3;
	transitionMatrix[1][0] = 0.4;
	transitionMatrix[1][1] = 0.6;

	double **emissionMatrix = new double*[2];
	for (int i = 0; i < 2; i++)
		emissionMatrix[i] = new double[3];
	emissionMatrix[0][0] = 0.5;
	emissionMatrix[0][1] = 0.4;
	emissionMatrix[0][2] = 0.1;
	emissionMatrix[1][0] = 0.1;
	emissionMatrix[1][1] = 0.3;
	emissionMatrix[1][2] = 0.6;

	double *initialProbabilities = new double[2];
	initialProbabilities[0]=0.6;
	initialProbabilities[1] = 0.4;

	HMM<string, string> hmm = HMM<string, string>(states, 2, observations, 3, transitionMatrix, emissionMatrix, initialProbabilities);

	string *testObservationSequence = new string[3];
	testObservationSequence[0]="normal";
	testObservationSequence[1]="cold";
	testObservationSequence[2]="dizzy";

	cout << "Testing for observation sequence: { " << testObservationSequence[0] << testObservationSequence[1] << testObservationSequence[2] << "}\n\n";

	cout << "Forward algorithm: " << hmm.forwardAlgorithm(testObservationSequence,3) << "\n";
	cout << "Backward algorthm: " << hmm.backwardAlgorithm(testObservationSequence,3) << "\n";

	string *result = new string[3];
	result = hmm.viterbiAlgorithm(testObservationSequence,3);
	
	cout << "Viterbi algorithm: { " << result[0] << result[1] << result[2] << "}\n";

	result = hmm.forwardBackwardAlgorithm(testObservationSequence, 3);

	cout << "ForwardBackward algorithm: { " << result[0] << result[1] << result[2] << "}\n";

	delete[] testObservationSequence;
	delete[] result;

	return 0;
}