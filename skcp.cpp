#include "SCPv.hpp"
#include "Random.hpp"
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>
using namespace std;


vector<int> SKCC;
vector<int> COST;
vector<int> SCORE;
vector<int> TIME;

const int max_iter = 1000;

// csに含まれない列から最大スコアのものを選んで返す
int get_column_maxscore(SCPinstance &inst,
                        SCPsolution& cs,
                        Rand& rnd)
{
  std::vector<int> maxCols;
  double maxScore = 0.0;
  int maxc = 0;
  double scw = 0.0;

  for (int c = 0; c < inst.numColumns; c++)
  {
    if (cs.SOLUTION[c]) { continue; }

    scw = (double)SCORE[c]/(double)inst.Weight[c];
    // 最大スコアの列をチェック
    if (maxScore < scw)
    {
      maxScore = SCORE[c];
      maxCols.clear();
      maxCols.push_back(c);
    }
    else if (maxScore == scw)
      maxCols.push_back(c);
  } // End for c

  if (maxCols.size() == 1) maxc = maxCols[0];
  else
  {
    int j = rnd(0, maxCols.size() - 1);
    maxc = maxCols[j];
  }

  return maxc;
}


int get_add_rule(SCPinstance &inst,
		 SCPsolution& cs,
		 Rand& rnd)
{
  std::vector<int> maxCols;
  double maxScore = 0.0;
  int retc = 0;
  double scw = 0.0;
  int oldest_time = max_iter + 1;

  for (int c = 0; c < inst.numColumns; c++)
  {
    if (cs.SOLUTION[c] == 1) { continue; }
    if (SKCC[c] == 0) { continue; }

    scw = (double)SCORE[c]/(double)inst.Weight[c];
    // 最大スコアの列をチェック
    if (maxScore < scw)
    {
      maxScore = SCORE[c];
      maxCols.clear();
      maxCols.push_back(c);
    }
    else if (maxScore == scw)
      maxCols.push_back(c);
  } // End for c

  if (maxCols.size() == 1) retc = maxCols[0];
  else
  {
    // とりあえずTIMEが同点のときは考えてない
    for (int c : maxCols) {
      if (TIME[c] < oldest_time) {
	oldest_time = TIME[c];
	retc = c;
      }
    }
    //int j = rnd(0, maxCols.size() - 1);
    //maxc = maxCols[j];
  }

  return retc;
} // add_rule


// REMOVE-RULE
int get_remove_rule(SCPinstance &inst,
		    SCPsolution& cs,
		    Rand& rnd)
{
  std::vector<int> maxCols;
  double maxScore = numeric_limits<int>::min();
  int retc = 0;
  double scw = 0.0;

  int oldest_time = max_iter + 1;

  for (int c : cs.CS) {
    scw = (double)SCORE[c]/(double)inst.Weight[c];

    // 最大スコアの列をチェック
    if (maxScore < scw)  {
      maxScore = SCORE[c];
      maxCols.clear();
      maxCols.push_back(c);
    }
    else if (maxScore == scw)
      maxCols.push_back(c);
  } // End for c

  if (maxCols.size() == 1) retc = maxCols[0];
  else
  {
    // とりあえずTIMEが同点のときは考えてない
    for (int c : maxCols) {
      if (TIME[c] < oldest_time) {
	oldest_time = TIME[c];
	retc = c;
      }
    }
    // int j = rnd(0, maxCols.size() - 1);
    // maxc = maxCols[j];
  }

  return retc;
}




// 配列の順序をランダムに入れ替える
void random_permutation(vector<int>& A, Rand& rnd)
{
  int j;
  int n = A.size();
  for (int i = 0; i < n-1; ++i)
  {
    j = rnd(i, n-1);
    swap(A[i], A[j]);
  }
}


void add_update_score(SCPinstance& inst, SCPsolution& cs, int c)
{
  SCORE[c] = -SCORE[c];

  for (int r : inst.ColEntries[c])  {
    // r行が初めてカバーされたら，rを含む行のスコアを減少
    if (cs.COVERED[r] == cs.K) {
      for (int rc : inst.RowCovers[r]) {
	if (rc != c) SCORE[rc] -= COST[r];
      } // End: for rc
    } // End if covered[r] == K

      // r行がK+1回カバーされたら，rを含むCSの要素のスコアを増加
    if (cs.COVERED[r] == cs.K+1) {
      for (int rc : inst.RowCovers[r]) {
	if (cs.SOLUTION[rc]) {
	  SCORE[rc] += COST[r];
	}
      }
    } // End if covered[r] == K+1
  } // end for r
}

void remove_update_score(SCPinstance& inst, SCPsolution& cs, int c)
{
  SCORE[c] = -SCORE[c];

  for (int r : inst.ColEntries[c])  {
    // r行がK回カバーされなくなったら，rを含む行のスコアを増加
    if (cs.COVERED[r] == cs.K-1) {
	for (int rc : inst.RowCovers[r]) {
	  if (rc != c && cs.SOLUTION[rc] == 0) SCORE[rc] += COST[r];
	} // End: for rc
      } // End if covered[r] == K-1

      // r行がK回カバーされたら，rを含むCSの列のスコアが減る
      if (cs.COVERED[r] == cs.K) {
	for (int rc : inst.RowCovers[r]) {
	  if (cs.SOLUTION[rc]) SCORE[rc] -= COST[r];
	}
      } // End if covered[r] == K
    } // end for r
} // end remove_update_score


// 列colの近傍のSKCCを1にする
void update_SKCC(SCPinstance& inst, int col)
{
  for (int c : inst.Neighborhood[col]) {
    SKCC[c] = 1;
  }
}



// 貪欲法：スコア最大の列をK列選ぶ
// 引数の cs に結果が入る
void greedy_construction(SCPinstance &inst,
                         SCPsolution &cs,
                         Rand &rnd)
{
  int maxc;

  cs.initialize(inst);
  while (cs.num_Cover < inst.numRows)
  {
    maxc = get_column_maxscore(inst, cs, rnd);
    cs.add_column(inst, maxc);

    // スコア更新
    add_update_score(inst, cs, maxc);
    // End update score
  } // End while num_Cover
}


SCPsolution DLL_com(SCPinstance& inst, int k, Rand& rnd)
{
  SCPsolution CS(inst, k);
  SCPsolution CSbest(inst, k);

  for (int j = 0; j < inst.numColumns; j++) {
    SCORE[j] = inst.ColEntries[j].size();
  }

  greedy_construction(inst, CS, rnd);
  CSbest = CS;

  int remove_col;

  for (int iter = 0; iter < max_iter; iter++) {
    cout << "Iter: " << iter;
    cout << " " << CSbest.totalWeight << " " << CS.totalWeight << " "
         << CS.num_Cover << " " << CS.CS.size() << " ";

    // 実行可能解が見つかったら更新
    if (CS.num_Cover == inst.numRows) {
      CSbest = CS;
      remove_col = get_remove_rule(inst, CS, rnd);
      CS.remove_column(inst, remove_col);
      TIME[remove_col] = iter;

      // update SCORE
      remove_update_score(inst, CS, remove_col);

      // update SKCC
      SKCC[remove_col] = 0;
      update_SKCC(inst, remove_col);
      // end update SKCC

      cout << " continue" << endl;
      continue;
    }

    // CS が実行可能でない場合
    // 1列削除する
    remove_col = get_remove_rule(inst, CS, rnd);
    cout << " Remove " << remove_col << " ";
    CS.remove_column(inst, remove_col);
    TIME[remove_col] = iter;
    remove_update_score(inst, CS, remove_col);

    // update SKCC
    SKCC[remove_col] = 0;
    update_SKCC(inst, remove_col);
    // end update SKCC


    int add_col;
    // 実行可能になるまで追加
    while (CS.num_Cover < inst.numRows) {
      add_col = get_add_rule(inst, CS, rnd);
      CS.add_column(inst, add_col);

      // 追加した結果が悪い解ならやめてやり直す
      if (CS.totalWeight >= CSbest.totalWeight) {
	CS.remove_column(inst, add_col);
	break;
      }

      add_update_score(inst, CS, add_col);
      // update SKCC
      update_SKCC(inst, add_col);
      // end update SKCC
      TIME[add_col] = iter;

      // update COST and SCORE
      for (int r = 0; r < inst.numRows; r++) {
	if (CS.COVERED[r] == 0) {
	  COST[r]++;
	  for (int rc : inst.RowCovers[r])  SCORE[rc]++;
	}
      }
    } // end while


    cout << endl;
  } // End iter

  return CSbest;
}


// メイン関数
int main(int argc, char** argv)
{
  //コマンドライン引数の数が少なければ強制終了
  if (argc < 3){
    cout << "Usage: ./command filename K(int)" << endl;
    return 0;
  }
  char *FileName = argv[1];
  FILE *SourceFile = fopen(FileName,"r");
  int K = atoi(argv[2]);

  // SCPのインスタンスを読み込む
  SCPinstance  instance(SourceFile);

  // initialize
  for (int i = 0; i < instance.numColumns; i++) {
    SKCC.push_back(1);
    SCORE.push_back(0);
    TIME.push_back(0);
  }
  for (int i = 0; i < instance.numRows; i++) {
    COST.push_back(1);
  }

  int seed = 0;
  Rand rnd;
  rnd.seed(seed);
  // End Initialize;


  SCPsolution CSbest(instance, K);

  CSbest = DLL_com(instance, K, rnd);

  CSbest.print_solution();
  cout << CSbest.totalWeight << " " << CSbest.num_Cover << endl;

  return 0;
}
