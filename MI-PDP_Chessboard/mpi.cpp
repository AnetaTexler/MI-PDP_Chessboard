#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "dirent.h"
#include <vector>
#include <utility>
#include <algorithm>
#include "math.h"
#include <chrono>
#include <omp.h>
#include "mpi.h"

#define MPI_PAR
#ifdef MPI_PAR

using namespace std;

class Task
{
public:
	int chessboardLen;
	int upperBound;
	vector<pair<int, int>> placedFiguresCoordinates;
	vector<pair<int, int>> resultMoves;
	pair<int, int> horseCoordinates;

	string toString()
	{
		string str = "---TASK------------------------------------------------------------\n";
		str += "Chessboard size: " + to_string(chessboardLen) + "x" + to_string(chessboardLen) + "\n";
		str += "Placed figures:  ";
		for (pair<int, int> coordinates : placedFiguresCoordinates)
			str += "(" + to_string(coordinates.first) + "," + to_string(coordinates.second) + ") ";
		str += "\n";
		str += "Horse position:  (" + to_string(horseCoordinates.first) + "," + to_string(horseCoordinates.second) + ")\n";
		str += "-------------------------------------------------------------------\n";
		return str;
	}
};

vector<pair<int, int>> g_resultMoves; // gives actual depth and number of moves (cost) 
vector<Task> g_taskPool;
vector<pair<int, int>> g_placedFiguresStart; // for preservation of start status

// Parse a file
void parseFile(ifstream & file, Task & inputData)//int & chessboardLen, int & upperBound, int & placedFiguresCnt, vector<pair<int, int>> & placedFiguresCoordinates, pair<int, int> & horseCoordinates)
{
	string line;
	getline(file, line);
	istringstream iss(line);
	iss >> inputData.chessboardLen >> inputData.upperBound;

	for (int r = 0; r < inputData.chessboardLen; r++)
	{
		getline(file, line);
		const char * x = line.c_str();
		for (int c = 0; c < inputData.chessboardLen; c++)
		{
			if (x[c] == '1')
				inputData.placedFiguresCoordinates.push_back(make_pair(r, c));
			if (x[c] == '3')
				inputData.horseCoordinates = make_pair(r, c);
		}
	}

	g_placedFiguresStart = inputData.placedFiguresCoordinates;
}

// Check if chessboard box contains a figure
bool containsFigure(const pair<int, int> & coordinatesToCheck, const vector<pair<int, int>> & placedFiguresCoordinates)
{
	if (find(placedFiguresCoordinates.begin(), placedFiguresCoordinates.end(), coordinatesToCheck) != placedFiguresCoordinates.end())
		return true;
	else
		return false;
}

// Return success rate of the specific move - on the box X: (SR(X) = 8*containsFigure(X) - distanceFromXToTheClosestFigure)
int getMoveSuccessRate(const pair<int, int> & moveCoordinates, const vector<pair<int, int>> & placedFiguresCoordinates, const int chessboardLen)
{
	int successRate = -(2 * chessboardLen - 2); // the biggest possible distance

	if (containsFigure(moveCoordinates, placedFiguresCoordinates))
		successRate = 8 * 1 - 0;
	else
	{
		for (pair<int, int> figCoordinates : placedFiguresCoordinates)
		{
			int tmp = 8 * 0 - (abs(moveCoordinates.first - figCoordinates.first) + abs(moveCoordinates.second - figCoordinates.second));
			if (tmp > successRate)
				successRate = tmp;
		}
	}

	return successRate;
}

// Try all possible moves if are valid (not out of chessboard size)
vector<pair<int, int>> getValidMoves(const int & chessboardLen, const pair<int, int> & horseCoordinates)
{
	vector<pair<int, int>> validMoves;

	if (horseCoordinates.first - 2 >= 0 && horseCoordinates.second - 1 >= 0)
		validMoves.push_back(make_pair(horseCoordinates.first - 2, horseCoordinates.second - 1));
	if (horseCoordinates.first - 1 >= 0 && horseCoordinates.second - 2 >= 0)
		validMoves.push_back(make_pair(horseCoordinates.first - 1, horseCoordinates.second - 2));
	if (horseCoordinates.first - 2 >= 0 && horseCoordinates.second + 1 < chessboardLen)
		validMoves.push_back(make_pair(horseCoordinates.first - 2, horseCoordinates.second + 1));
	if (horseCoordinates.first + 1 < chessboardLen && horseCoordinates.second - 2 >= 0)
		validMoves.push_back(make_pair(horseCoordinates.first + 1, horseCoordinates.second - 2));
	if (horseCoordinates.first - 1 >= 0 && horseCoordinates.second + 2 < chessboardLen)
		validMoves.push_back(make_pair(horseCoordinates.first - 1, horseCoordinates.second + 2));
	if (horseCoordinates.first + 2 < chessboardLen && horseCoordinates.second - 1 >= 0)
		validMoves.push_back(make_pair(horseCoordinates.first + 2, horseCoordinates.second - 1));
	if (horseCoordinates.first + 1 < chessboardLen && horseCoordinates.second + 2 < chessboardLen)
		validMoves.push_back(make_pair(horseCoordinates.first + 1, horseCoordinates.second + 2));
	if (horseCoordinates.first + 2 < chessboardLen && horseCoordinates.second + 1 < chessboardLen)
		validMoves.push_back(make_pair(horseCoordinates.first + 2, horseCoordinates.second + 1));

	return validMoves;
}

// Driver function to sort the vector elements by first element of pair in descending order
bool sortDesc(const pair<int, pair<int, int>> &a, const pair<int, pair<int, int>> &b)
{
	return (a.first > b.first);
}

// Sequential recursive function for search in state space
void solveInstanceSEQ(const int chessboardLen, const int upperBound, vector<pair<int, int>> placedFiguresCoordinates,
	pair<int, int> horseCoordinates, vector<pair<int, int>> resultMoves)
{

		if (placedFiguresCoordinates.empty()) // all figures removed by horse
		{
			if (g_resultMoves.size() == 0)
				g_resultMoves = resultMoves;
			else
				g_resultMoves = (resultMoves.size() < g_resultMoves.size()) ? resultMoves : g_resultMoves;
			return;
		}

		// needless to search in the same or higher depth than actual minimum ( = best cost = g_resultMoves.size()) OR upperBound is achieved
		if ((g_resultMoves.size() != 0 && (resultMoves.size() + placedFiguresCoordinates.size() >= g_resultMoves.size())) || resultMoves.size() > upperBound)
			return;

	vector<pair<int, int>> validMoves = getValidMoves(chessboardLen, horseCoordinates);
	vector<pair<int, pair<int, int>>> ratedValidMoves;

	// count move success rate for each valid move
	for (pair<int, int> validMove : validMoves)
	{
		int successRate = getMoveSuccessRate(validMove, placedFiguresCoordinates, chessboardLen);
		ratedValidMoves.push_back(make_pair(successRate, validMove));
	}

	sort(ratedValidMoves.begin(), ratedValidMoves.end(), sortDesc); // the closest figures first

	if (resultMoves.size() == 0)
		resultMoves.push_back(horseCoordinates); // start horse position

	for (pair<int, pair<int, int>> ratedValidMove : ratedValidMoves)
	{
		vector<pair<int, int>> placedFiguresCoordinatesCopy = placedFiguresCoordinates;
		vector<pair<int, int>> resultMovesCopy = resultMoves;

		// remove figure from vector placedFiguresCoordinates
		if (ratedValidMove.first == 8)
			placedFiguresCoordinatesCopy.erase(remove(placedFiguresCoordinatesCopy.begin(), placedFiguresCoordinatesCopy.end(), ratedValidMove.second), placedFiguresCoordinatesCopy.end());

		// add move if is not the same like the penultimate move (that would mean the loop until the upper bound is achieved)
		resultMovesCopy.push_back(ratedValidMove.second);

		// RECURSION CALL
		solveInstanceSEQ(chessboardLen, upperBound, placedFiguresCoordinatesCopy, ratedValidMove.second, resultMovesCopy);
	}
}

bool taskExistsInPool(Task taskToCheck)
{
	// 2 tasks with the same position of horse and with 0 or the same picked figures are duplicated
	for (Task existingTask : g_taskPool)
	{
		if (existingTask.horseCoordinates == taskToCheck.horseCoordinates) // same horse position
		{
			if (g_placedFiguresStart == taskToCheck.placedFiguresCoordinates || existingTask.placedFiguresCoordinates == taskToCheck.placedFiguresCoordinates) // 0 picked OR same figures picked
				return true;
		}
	}

	return false;
}

void GenerateTaskPool(Task inputTask)
{
	if (inputTask.resultMoves.size() == 5) // tree depth 
	{
		if (!taskExistsInPool(inputTask)) // reduce duplicated tasks
			g_taskPool.push_back(inputTask);
		return;
	}

	vector<pair<int, int>> validMoves = getValidMoves(inputTask.chessboardLen, inputTask.horseCoordinates);

	if (inputTask.resultMoves.size() == 0)
		inputTask.resultMoves.push_back(inputTask.horseCoordinates); // start horse position

	for (pair<int, int> validMove : validMoves)
	{
		Task taskCopy = inputTask;
		taskCopy.resultMoves.push_back(validMove); // add move
		taskCopy.horseCoordinates = validMove; // new position of horse
		if (containsFigure(validMove, taskCopy.placedFiguresCoordinates)) // if position contains a figure, remove it
			taskCopy.placedFiguresCoordinates.erase(remove(taskCopy.placedFiguresCoordinates.begin(), taskCopy.placedFiguresCoordinates.end(), validMove),
				taskCopy.placedFiguresCoordinates.end());

		GenerateTaskPool(taskCopy);
	}
}

vector<int> serialize(vector<pair<int, int>> & input)
{
	vector<int> result;
	for (pair<int, int> pair : input)
	{
		result.push_back(pair.first);
		result.push_back(pair.second);
	}
	return result;
}

vector<pair<int, int>> deserialize(int * input, int size)
{
	vector<pair<int, int>> result;
	for (int i = 0; i < size - 1; i += 2)
	{
		result.push_back(make_pair(input[i], input[i + 1]));
	}
	return result;
}



int main(int argc, char *argv[])
{
	string dirPath = "inputFiles";
	DIR * dir;
	struct dirent * entry;
	Task inputTask;

	dir = opendir(dirPath.c_str());
	if (dir == NULL) return 1;
	
	while (entry = readdir(dir))
	{
		if ((string)(entry->d_name) == "." || (string)(entry->d_name) == "..") continue;
		ifstream ifs(dirPath + "/" + (string)(entry->d_name));

		parseFile(ifs, inputTask);
		auto start = chrono::system_clock::now();

		// MASTER - SLAVE parallelism --------------------------

		// int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
			// buf - ukazatel na posilana data
			// count - pocet posilanych polozek
			// datatype - datovy typ posilanych dat
			// dest - cislo ciloveho procesu
			// tag - znacka zpravy
			// comm - platny MPI komunikator (MPI_COMM_WORLD)

		// int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
			// buf - ukazatel na buffer pro prijmana data
			// count - maximaln pocet prijmanych polozek
			// datatype - urcuje datovy typ prijmanych dat
			// source - cislo zdrojoveho procesu
			// tag - znacka zpravy
			// comm - platny MPI komunikator (MPI COMM WORLD)
			// status - ukazatel na stavovy objekt
			
		const int MPI_MASTER_ID = 0;
		const int MPI_TAG_WORK_DEMAND = 1;
		const int MPI_TAG_WORK_COMMAND = 2;
		const int MPI_TAG_RESULT = 3;
		const int MPI_TAG_NO_RESULT = 4;
		const int MPI_TAG_END = 5;

		int processID;
		int processCount;

		MPI_Init(NULL, NULL);
		MPI_Comm_rank(MPI_COMM_WORLD, &processID);
		MPI_Comm_size(MPI_COMM_WORLD, &processCount);

		GenerateTaskPool(inputTask); // hack
		//cout << "TaskPoolGenerated: " << processID << endl;
		//MPI_Barrier(MPI_COMM_WORLD);

		// define master process - e.g. process with id = 0
		
		// MASTER (0)
		if (processID == 0)
		{
			// GenerateTaskPool(inputTask);

			int currentPoolPosition = 0;

			while (true)
			{
				int flag; // message arrival flag
				MPI_Status status;
				MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_WORK_DEMAND, MPI_COMM_WORLD, &flag, &status); // message arrival test without accepting it (non-block version)

				if (flag) // message arrived - master assigns work
				{
					int message;
					MPI_Status status2;
					MPI_Recv(&message, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status2);
					int slaveID = status2.MPI_SOURCE;

					if (currentPoolPosition == g_taskPool.size() - 1) // remains the last task to process
					{
						MPI_Send(&currentPoolPosition, 1, MPI_INT, slaveID, MPI_TAG_WORK_COMMAND, MPI_COMM_WORLD); // send to correct slave! -> slaveID
						cout << "MASTER - sends to SLAVE " << slaveID << ": task pool position: " << currentPoolPosition << ", the last task" << endl;
						break; // go to collect results
					}

					MPI_Send(&currentPoolPosition, 1, MPI_INT, slaveID, MPI_TAG_WORK_COMMAND, MPI_COMM_WORLD);
					cout << "MASTER - sends to SLAVE " << slaveID << ": task pool position: " << currentPoolPosition << endl;
					
					currentPoolPosition++;
				}
				else // no message arrived - master works
				{
					//continue;
					if (currentPoolPosition == g_taskPool.size() - 1) // remains the last task to process
					{
						cout << "MASTER - solves the last task:" << endl << g_taskPool[currentPoolPosition].toString() << endl;
						solveInstanceSEQ(g_taskPool[currentPoolPosition].chessboardLen,
										 g_taskPool[currentPoolPosition].upperBound,
										 g_taskPool[currentPoolPosition].placedFiguresCoordinates,
										 g_taskPool[currentPoolPosition].horseCoordinates,
										 g_taskPool[currentPoolPosition].resultMoves);
						break; // go to collect results
					}
					cout << "MASTER - solves: pool position: " << currentPoolPosition << "/" << g_taskPool.size() - 1 << endl << g_taskPool[currentPoolPosition].toString() << endl;
					solveInstanceSEQ(g_taskPool[currentPoolPosition].chessboardLen,
									 g_taskPool[currentPoolPosition].upperBound,
									 g_taskPool[currentPoolPosition].placedFiguresCoordinates,
									 g_taskPool[currentPoolPosition].horseCoordinates,
									 g_taskPool[currentPoolPosition].resultMoves);

					currentPoolPosition++;
				}
			}
		}

		// SLAVES (from 1 to processCount - 1)
		else
		{
			while (true) 
			{
				//cout << "SLAVE " << processID << " MPI_TAG_WORK_DEMAND - will send" << endl;
				int currentPoolPosition;
				MPI_Status status;
				MPI_Send(0, 0, MPI_INT, MPI_MASTER_ID, MPI_TAG_WORK_DEMAND, MPI_COMM_WORLD); // send request for work to master
				//cout << "SLAVE " << processID << " MPI_TAG_WORK_DEMAND - was sent" << endl;
				MPI_Recv(&currentPoolPosition, 1, MPI_INT, MPI_MASTER_ID, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // receive any answer from master
				//cout << "SLAVE " << processID << " MPI_TAG_WORK_DEMAND - response: " << status.MPI_TAG << endl;

				if (status.MPI_TAG == MPI_TAG_WORK_COMMAND) // work assigned by master
				{
					cout << "SLAVE " << processID << " - command to solve: pool position: " << currentPoolPosition << "/" << g_taskPool.size() - 1 << endl << g_taskPool[currentPoolPosition].toString() << endl;
					solveInstanceSEQ(g_taskPool[currentPoolPosition].chessboardLen, 
									 g_taskPool[currentPoolPosition].upperBound, 
									 g_taskPool[currentPoolPosition].placedFiguresCoordinates, 
									 g_taskPool[currentPoolPosition].horseCoordinates,
									 g_taskPool[currentPoolPosition].resultMoves);
					//cout << "SLAVE " << processID << " - solved: pool position: " << currentPoolPosition << "/" << g_taskPool.size() - 1 << endl << g_taskPool[currentPoolPosition].toString() << endl;
				}
				else if (status.MPI_TAG == MPI_TAG_END) // no more work, master wants result from slave
				{
					if (g_resultMoves.empty()) // no result
					{
						MPI_Send(0, 0, MPI_INT, MPI_MASTER_ID, MPI_TAG_NO_RESULT, MPI_COMM_WORLD); 
					}
					else
					{
						vector<int> serializedResultMoves = serialize(g_resultMoves);
						int tmpSize = serializedResultMoves.size();
						MPI_Send(&tmpSize, 1, MPI_INT, MPI_MASTER_ID, MPI_TAG_RESULT, MPI_COMM_WORLD); // size
						MPI_Send(serializedResultMoves.data(), serializedResultMoves.size(), MPI_INT, MPI_MASTER_ID, MPI_TAG_RESULT, MPI_COMM_WORLD); // array
					}
					break;
				}
				else
				{
					cout << "ERROR: SLAVE " << processID << " - unexpected MPI_TAG: " << status.MPI_TAG << endl;
					abort();
				}
			}
		}

		// MASTER - collects results from slaves
		if (processID == 0) 
		{
			int collectedResultsCnt = 1; // 1 for MASTER
			for (collectedResultsCnt = 1; collectedResultsCnt < processCount;)
			{
				int size;
				MPI_Status status;
				MPI_Recv(&size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				
				if (status.MPI_TAG == MPI_TAG_WORK_DEMAND) // send END_TAG to slaves who want work
				{
					cout << "MASTER - sends END_TAG to SLAVE: " << status.MPI_SOURCE << endl;
					MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, MPI_TAG_END, MPI_COMM_WORLD);
					continue;
				}
				else if (status.MPI_TAG == MPI_TAG_RESULT) // slave send result, master checks if the count of moves is lower than in current best result 
				{
					int * serializedResultMoves = new int[size];
					MPI_Recv(serializedResultMoves, size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

					vector<pair<int, int>> resultMoves = deserialize(serializedResultMoves, size);

					if (g_resultMoves.empty() || resultMoves.size() < g_resultMoves.size())
					{
						g_resultMoves = resultMoves;
					}

					cout << "SLAVE " << status.MPI_SOURCE << " returned result size: " << resultMoves.size() << endl;
					collectedResultsCnt++;
				}
				else if (status.MPI_TAG == MPI_TAG_NO_RESULT)
				{
					cout << "SLAVE " << status.MPI_SOURCE << " returned no result." << endl;
					collectedResultsCnt++;
				}
				else
				{
					cout << "ERROR: MASTER - unexpected MPI_TAG: " << status.MPI_TAG << endl;
					abort();
				}
			}

			auto end = chrono::system_clock::now();
			long long millis = chrono::duration_cast<chrono::milliseconds>(end - start).count();
			cout << "=====================================================================================" << endl;
			cout << "SOLUTION FOR " << (string)(entry->d_name) << ":" << endl;
			cout << g_taskPool.size() << " tasks" << endl;

			if (g_resultMoves.size() == 0)
				cout << "Not found." << endl << endl;
			else
			{
				cout << "Time: " << millis / 1000.0 << "s, Count: " << g_resultMoves.size() - 1 << ", Moves: ";
				for (pair<int, int> resultMove : g_resultMoves)
				{
					cout << "(" << resultMove.first << "," << resultMove.second << ")";
					if (containsFigure(resultMove, inputTask.placedFiguresCoordinates))
						cout << "* ";
					else cout << " ";
				}
			}
			cout << endl << "=====================================================================================" << endl << endl;
		}

		if (processID)
			cout << "SLAVE " << processID << " - finalize." << endl;
		else
			cout << "MASTER - finalize.";
		MPI_Finalize();
		
		g_resultMoves.clear();
		g_taskPool.clear();
		g_placedFiguresStart.clear();
		inputTask.placedFiguresCoordinates.clear();

		ifs.close();
	}

	closedir(dir);

	cin.get();

	return 0;
}

#endif
