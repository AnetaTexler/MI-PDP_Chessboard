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

using namespace std;

#define OMP_TASK_PAR
#ifdef OMP_TASK_PAR

// COMPILE COMMAND
// g++ -fopenmp main.cpp 

vector<pair<int, int>> g_resultMoves; // gives actual depth and number of moves (cost) 

// Parse a file
void parseFile(ifstream & file, int & chessboardLen, int & upperBound, int & placedFiguresCnt, vector<pair<int, int>> & placedFiguresCoordinates,
			   pair<int, int> & horseCoordinates)
{
	string line;
	getline(file, line);
	istringstream iss(line);
	iss >> chessboardLen >> upperBound;

	for (int r = 0; r < chessboardLen; r++)
	{
		getline(file, line);
		const char * x = line.c_str();
		for (int c = 0; c < chessboardLen; c++)
		{
			if (x[c] == '1')
			{
				placedFiguresCoordinates.push_back(make_pair(r, c));
				placedFiguresCnt++;
			}
			if (x[c] == '3')
				horseCoordinates = make_pair(r, c);
		}
	}

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
void solveInstanceSEQ(const int chessboardLen, const int upperBound, const int placedFiguresCnt, vector<pair<int, int>> placedFiguresCoordinates,
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

	// needless to search in the same or higher depth than actual minimum ( = best cost = RESULT_MOVES.size()) OR upperBound is achieved
	if ((g_resultMoves.size() != 0 && (resultMoves.size() + placedFiguresCoordinates.size() >= g_resultMoves.size()))
		|| resultMoves.size() > upperBound)
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
		solveInstanceSEQ(chessboardLen, upperBound, placedFiguresCnt, placedFiguresCoordinatesCopy, ratedValidMove.second, resultMovesCopy);
	}
}

// Parallel recursive function for search in state space
void solveInstanceOMP(const int chessboardLen, const int upperBound, const int placedFiguresCnt, vector<pair<int, int>> placedFiguresCoordinates,
					  pair<int, int> horseCoordinates, vector<pair<int, int>> resultMoves)
{
	bool exit = false;
	#pragma omp critical
	{
		if (placedFiguresCoordinates.empty()) // all figures removed by horse
		{
			if (g_resultMoves.size() == 0)
				g_resultMoves = resultMoves;
			else
				g_resultMoves = (resultMoves.size() < g_resultMoves.size()) ? resultMoves : g_resultMoves;
			//return;
			exit = true;
		}

		// needless to search in the same or higher depth than actual minimum ( = best cost = g_resultMoves.size()) OR upperBound is achieved
		if ((g_resultMoves.size() != 0 && (resultMoves.size() + placedFiguresCoordinates.size() >= g_resultMoves.size())) || resultMoves.size() > upperBound)
			//return;
			exit = true; 
	}

	if (exit)
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
		if (resultMovesCopy.size() < 5)
		{
			#pragma omp task
			solveInstanceOMP(chessboardLen, upperBound, placedFiguresCnt, placedFiguresCoordinatesCopy, ratedValidMove.second, resultMovesCopy);
		}
		else
			solveInstanceSEQ(chessboardLen, upperBound, placedFiguresCnt, placedFiguresCoordinatesCopy, ratedValidMove.second, resultMovesCopy);
	}
}

// Parallel start function
void openMPSolveInstance(const int chessboardLen, const int upperBound, const int placedFiguresCnt, vector<pair<int, int>> placedFiguresCoordinates,
						 pair<int, int> horseCoordinates, vector<pair<int, int>> resultMoves)
{
	#pragma omp parallel
	{
		#pragma omp single
		solveInstanceOMP(chessboardLen, upperBound, placedFiguresCnt, placedFiguresCoordinates, horseCoordinates, resultMoves);
	}
}


int main()
{
	string dirPath = "inputFiles";
	DIR * dir;
	struct dirent * entry;
	int chessboardLen = 0, upperBound = 0, placedFiguresCnt = 0;
	vector<pair<int, int>> placedFiguresCoordinates, resultMoves;
	pair<int, int> horseCoordinates(0, 0);

	dir = opendir(dirPath.c_str());
	if (dir != NULL)
	{
		while (entry = readdir(dir))
		{
			if ((string)(entry->d_name) == "." || (string)(entry->d_name) == "..") continue;
			ifstream ifs(dirPath + "\\" + (string)(entry->d_name));

			parseFile(ifs, chessboardLen, upperBound, placedFiguresCnt, placedFiguresCoordinates, horseCoordinates);

			auto start = chrono::system_clock::now();

			openMPSolveInstance(chessboardLen, upperBound, placedFiguresCnt, placedFiguresCoordinates, horseCoordinates, resultMoves);
			//solveInstanceSEQ(chessboardLen, upperBound, placedFiguresCnt, placedFiguresCoordinates, horseCoordinates, resultMoves);

			auto end = chrono::system_clock::now();
			long long millis = chrono::duration_cast<chrono::milliseconds>(end - start).count();

			cout << "SOLUTION FOR " << (string)(entry->d_name) << ":" << endl;
			if (g_resultMoves.size() == 0)
				cout << "Not found." << endl << endl;
			else
			{
				cout << "Time: " << millis / 1000.0 << "s, Count: " << g_resultMoves.size() - 1 << ", Moves: ";
				for (pair<int, int> resultMove : g_resultMoves)
				{
					cout << "(" << resultMove.first << "," << resultMove.second << ")";
					if (containsFigure(resultMove, placedFiguresCoordinates))
						cout << "* ";
					else cout << " ";
				}
				cout << endl << endl;
			}

			g_resultMoves.clear();
			placedFiguresCnt = 0;
			placedFiguresCoordinates.clear();

			ifs.close();
		}

	}
	closedir(dir);

	cin.get();

	return 0;
}

#endif
