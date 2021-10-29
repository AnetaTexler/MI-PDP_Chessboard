# MI-PDP_Chessboard
## Parallel and distributed programming course

### Task: Knight on a chessboard

**Input data:**
* k = a natural number, k > 5, size of a chessbord C (kxk) 
* q = a natural number, q < k^2-1, number of pieces on a chessboard C 
* F[1..q] = an array of coordinates of pieces 
* J = a coordinate of a knight 

**Rules and game goal:**

At the beginning, there are q pieces and 1 knight (horse) on the square chessboard C. We will call this arrangement the initial configuration. One move is the knight's jump according to the chess rules. The aim of the game is to capture all pieces with a knight in a minimum number of moves so that only the knight himself remains on the board. 

**Output of the algorithm:**

The minimum sequence of knight's moves leading to the goal (i.e. state where the knight remains alone on the board). The output will be a list of coordinates where the knight moves with an asterisk mark of the squares where a piece was captured. 

**Definition:**

Function val(X) evaluates any square X of the chessboard C as follows: if the square X contains a piece, the function val(X) returns the value 1, otherwise it returns the value 0. The success of the knight's move to the square X can be defined as a function: 

U(X) = 8*val(X) - distance_from_X_to_nearest_piece

**Sequential algorithm:**

A BB-DFS type with the depth of the searched space limited by the upper bound (see below). There is always a solution.
The permissible end state is a situation where the knight himself remains on the board. 
The cost we minimize is the number of moves to get to the final state from the initial configuration.
The algorithm ends when the cost is equal to the lower bound (see below) or when the entire state space is searched to a depth given by the upper bound. 

Recommendation for implementation: when implementing a sequential algorithm, make successors X of the current square of knight Y in the order of decreasing success function U (X) of the move to square X. This actually jumps to the squares where you can immediately capture a piece and there is a change to capture a piece in minimum moves.
The lower bound on the cost of the solution is q, because the knight can pick up a maximum of one piece in one move, and therefore we cannot find a solution faster than in q moves. 

The upper bound on the cost of the solution is k * k-1, as this is the minimum number of moves when the knight visits all squares of the chessboard C. This bound is too pessimistic in some situations, so the recommended upper bound value will be given directly in the initial configuration data file.

Pruning tree branches: It does not make sense to look for a solution at the same or greater depth than the minimum already found (curr_min). Therefore: 

curr_depth + (q - num_of_captured_pieces) < curr_min,

else the branch is terminated.

**Parallel algorithm:** task and data parallelism (OpenMP), Master-Slave (MPI).
