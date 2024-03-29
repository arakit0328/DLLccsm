#include "SCPv.hpp"
#include "Random.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>
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
vector<int> TIMES;


int compute_score(SCPinstance& inst,
                  SCPsolution& cs,
                  int c)
{
  int sc = 0;
  for (int r : inst.ColEntries[c]) {
    if (cs.SOLUTION[c] && (cs.COVERED[r] == cs.K)) sc -= COST[r];
    else if (!cs.SOLUTION[c] && cs.COVERED[r] < cs.K) sc += COST[r];
  }
  return sc;
}


// csに含まれない列から最大スコアのものを選んで返す
int get_column_maxscore(SCPinstance &inst,
                        SCPsolution& cs,
                        vector<int>& score,
                        Rand& rnd)
{
  std::vector<int> maxCols;
  double maxScore = 0.0;
  int maxc = 0;
  double scw = 0.0;

  for (int c = 0; c < inst.numColumns; c++)
  {
    if (cs.SOLUTION[c]) { continue; }

    scw = (double)score[c]/(double)inst.Weight[c];
    // 最大スコアの列をチェック
    if (maxScore < scw)
    {
      maxScore = scw;
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
  int oldest_time = numeric_limits<int>::max();

  for (int c = 0; c < inst.numColumns; c++)
  {
    if (cs.SOLUTION[c] == 1) { continue; }
    if (SKCC[c] == 0) { continue; }


    scw = (double)SCORE[c]/(double)inst.Weight[c];
    // 最大スコアの列をチェック
    if (maxScore < scw)
    {
      maxScore = scw;
      maxCols.clear();
      maxCols.push_back(c);
    }
    else if (maxScore == scw)
      maxCols.push_back(c);
  } // End for c

  if (maxCols.size() == 1) retc = maxCols[0];
  else
  {
    for (int c : maxCols) {
      if (TIMES[c] < oldest_time) {
        oldest_time = TIMES[c];
        retc = c;
      }
    }
  }

  return retc;
} // add_rule


// REMOVE-RULE
int get_remove_rule(SCPinstance &inst,
		    SCPsolution& cs,
                    int iter,
		    Rand& rnd)
{
  std::vector<int> maxCols;
  double maxScore = numeric_limits<int>::min();
  int retc = 0;
  double scw = 0.0;

  int oldest_time = numeric_limits<int>::max();

  if (rnd() % 100 < 95)
  {
    for (int c : cs.CS) {
      if (TIMES[c] > 0 && TIMES[c] == iter - 1) continue;

      // Araki
      // スコアが0の列の取り扱い
      // すべての行をk回カバーしている場合のみ取り除く
      bool flg = false;
      if (SCORE[c] == 0) {
        for (int r : inst.ColEntries[c]) {
          if (cs.COVERED[r] < cs.K) {
            flg = true;
            break;
          }
        }
        if (flg) continue;
      }

      scw = (double)SCORE[c]/(double)inst.Weight[c];

      // 最大スコアの列をチェック
      if (maxScore < scw)  {
        maxScore = scw;
        maxCols.clear();
        maxCols.push_back(c);
      }
      else if (maxScore == scw)
        maxCols.push_back(c);
    } // End for c

    if (maxCols.size() == 1) retc = maxCols[0];
    else
    {
      for (int c : maxCols) {
        if (TIMES[c] < oldest_time) {
          oldest_time = TIMES[c];
          retc = c;
        }
      }
    }
  }
  else
  {
    // 5%
    int maxw = 0;
    for (int c : cs.CS) {
      if (TIMES[c] < oldest_time) {
        oldest_time = TIMES[c];
        maxw = inst.Weight[c];
        retc = c;
      }
      else if (TIMES[c] == oldest_time) {
        if (maxw < inst.Weight[c]) {
          maxw = inst.Weight[c];
          retc = c;
        }
      }
    }
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
  SCORE[c] = 0;
  for (int r : inst.ColEntries[c])
  {
    if (cs.COVERED[r] == cs.K) SCORE[c] -= COST[r];
  }

  for (int r : inst.ColEntries[c])
  {
    if (cs.COVERED[r] == cs.K)
    {
      for (int rc : inst.RowCovers[r]) {
        if (rc != c) SCORE[rc] -= COST[r];
      }
    }
    else if (cs.COVERED[r] == cs.K + 1)
    {
      for (int rc : inst.RowCovers[r])
      {
	if (cs.SOLUTION[rc] && rc != c) {
	  SCORE[rc] += COST[r];
	}
      }
    } // End if covered[r] == K+1
  } // end for r
}


void remove_update_score(SCPinstance& inst, SCPsolution& cs, int c)
{
  SCORE[c] = 0;
  for (int r : inst.ColEntries[c])
  {
    if (cs.COVERED[r] < cs.K) SCORE[c] += COST[r];
  }

  for (int r : inst.ColEntries[c])
  {
    // r行がK回カバーされなくなったら，rを含む行のスコアを増加
    if (cs.COVERED[r] == cs.K-1) {
	for (int rc : inst.RowCovers[r])
        {
	  if (rc != c)
          {
            SCORE[rc] += COST[r];
          }
	} // End: for rc
      }
    else if (cs.COVERED[r] == cs.K)
    {
      for (int rc : inst.RowCovers[r])
      {
        if (cs.SOLUTION[rc] && rc != c)
        {
          SCORE[rc] -= COST[r];
        }
      }
    }
  } // end for r
} // end remove_update_score


// 列colの近傍のSKCCを1にする
void update_SKCC(SCPinstance& inst, SCPsolution& cs, int col)
{
  for (int c : inst.Neighborhood[col]) {
      SKCC[c] = 1;
  }
  // for (int r : inst.ColEntries[col]) {
  //   if (cs.COVERED[r] < cs.K) {
  //     for (int rc : inst.RowCovers[r]) {
  //       if (rc != col) SKCC[rc] = 1;
  //     }
  //   }
  // }
}



// 貪欲法：スコア最大の列をK列選ぶ
// 引数の cs に結果が入る
SCPsolution greedy_construction(SCPinstance &inst,
                                int k,
                                Rand &rnd)
{
  vector<int> score(inst.numColumns, 0);

  for (int c = 0; c < inst.numColumns; c++)
    score[c] = inst.ColEntries[c].size();

  int ca;
  SCPsolution cs(inst, k);

  while (cs.num_Cover < inst.numRows)
  {
    ca = get_column_maxscore(inst, cs, score, rnd);
    cs.add_column(inst, ca);

    // スコア更新
    score[ca] = 0;
    for (int r : inst.ColEntries[ca])
    {
      if (cs.COVERED[r] == cs.K)
      {
        for (int rc : inst.RowCovers[r])
          if (!cs.SOLUTION[rc] && rc != ca) SCORE[rc]--;
      }
    } // end for r
  } // End while num_Cover

  return cs;
}




SCPsolution DLL_com(SCPinstance& inst, int k, int max_iter, Rand& rnd)
{
  SCPsolution CS(inst, k);
  SCPsolution CSbest(inst, k);

  vector<int> Freq(inst.numColumns, 0);

  CS = greedy_construction(inst, k, rnd);
  CSbest = CS;

  for (int c : CS.CS) TIMES[c] = 1;

  for (int c = 0; c < inst.numColumns; c++)
  {
    SCORE[c] = 0;
    if (CS.SOLUTION[c])
    {
      for (int r : inst.ColEntries[c])
        if (CS.COVERED[r] == k) SCORE[c] -= COST[r];
    }
  }

  int remove_col;

  for (int iter = 1; iter <= max_iter; iter++)
  {
    // cout << "Iter: " << iter;
    // cout << " " << CSbest.totalWeight << " " << CS.totalWeight << " " << CS.num_Cover << " " << CS.CS.size() << " ";



    // 実行可能解が見つかったら更新
    if (CS.num_Cover == inst.numRows) {
      CSbest = CS;
      remove_col = get_remove_rule(inst, CS, 0, rnd);

      // cout << " Remove " << remove_col << "(" << (double)SCORE[remove_col] / inst.Weight[remove_col] << ") ";

      CS.remove_column(inst, remove_col);
      TIMES[remove_col] = iter;
      //Freq[remove_col]++;

      // update SCORE
      remove_update_score(inst, CS, remove_col);

      // update SKCC
      SKCC[remove_col] = 0;
      update_SKCC(inst, CS, remove_col);
      // end update SKCC

      // cout << " continue" << endl;
      continue;
    }

    // CS が実行可能でない場合
    // 1列削除する
    remove_col = get_remove_rule(inst, CS, iter, rnd);
    // cout << " Remove " << remove_col << "(" << (double)SCORE[remove_col] / inst.Weight[remove_col] << ") ";
    CS.remove_column(inst, remove_col);
    TIMES[remove_col] = iter;
    remove_update_score(inst, CS, remove_col);
    //Freq[remove_col]++;

    // update SKCC
    SKCC[remove_col] = 0;
    update_SKCC(inst, CS, remove_col);
    // end update SKCC

    int add_col;
    //bool brk_flag = false;

    // 実行可能になるまで追加
    while (CS.num_Cover < inst.numRows) {
      add_col = get_add_rule(inst, CS, rnd);

      if (CS.totalWeight + inst.Weight[add_col] >= CSbest.totalWeight)
      {
        // 追加した結果が悪い解ならやめてやり直す
        //brk_flag = true;
        // cout << "Add " << add_col << " Break ";
	break;
      }
      else
      {
        // cout << "Add " << add_col << "(" << (double)SCORE[add_col] / inst.Weight[add_col] << ") ";
        CS.add_column(inst, add_col);
        add_update_score(inst, CS, add_col);

        update_SKCC(inst, CS, add_col);
        SKCC[add_col] = 0;

        TIMES[add_col] = iter;
        Freq[add_col]++;

        // Araki: COST reset
        // add_colを追加してK回カバーされた列のcostを1に戻してスコア再計算
        for (int r : inst.ColEntries[add_col])
        {
          if (COST[r] > max_iter / 10 && CS.COVERED[r] == CS.K)
          {
            COST[r] = 1;
            for (int rc : inst.RowCovers[r])
            {
              SCORE[rc] = compute_score(inst, CS, rc);
            }
          }
        }
      }

      // update COST and SCORE
      for (int r = 0; r < inst.numRows; r++)
      {
        if (CS.COVERED[r] < CS.K)
        {
          COST[r]++;
          for (int rc : inst.RowCovers[r])
          {
            if (!CS.SOLUTION[rc]) SCORE[rc]++;
          }
        }
      }

      // Check
      // for (int c = 0; c < inst.numColumns; c++) {
      //   int sc = compute_score(inst, CS, c);
      //   if (sc != SCORE[c]) {
      //     printf("Col %d (%d) sc = %d, but SCORE[%d] = %d\n", c, CS.SOLUTION[c], sc, c, SCORE[c]);
      //     exit(1);
      //   }
      // }

    } // end while CS.num_Cover

    // cout << endl;
  } // End iter


  // for (int c= 0; c < inst.numColumns; c++) {
  //   cout << c << " ";
  //   if (CSbest.SOLUTION[c]) cout << "+ ";
  //   else cout << "  ";
  //   cout << Freq[c] << endl;
  // }

  return CSbest;
}



bool check_solution(SCPinstance& inst, SCPsolution& cs)
{
  vector<int> cov(inst.numRows, 0);
  int tw = 0;

  for (int c : cs.CS) {
    tw += inst.Weight[c];
    for (int r : inst.ColEntries[c]) {
      cov[r]++;
    }
  }

  for (int r = 0; r < inst.numRows; r++) {
    if (cov[r] < cs.K) {
      cout << "This is not a feasible solution." << endl;
      return false;
    }
  }

  if (tw != cs.totalWeight) {
    cout << "Wrong totalWeight." << endl;
    return false;
  }

  return true;
}


// メイン関数
int main(int argc, char** argv)
{
  //コマンドライン引数の数が少なければ強制終了
  // if (argc < 3){
  //   cout << "Usage: ./command filename K(int)" << endl;
  //   return 0;
  // }
  // char *FileName = argv[1];
  // FILE *SourceFile = fopen(FileName,"r");
  // int K = atoi(argv[2]);

  //コマンドライン引数の数が少なければ強制終了
  if (argc < 2) {
    cout << "Usage: ./command filename" << endl;
    return 0;
  }

  char *FileName = argv[1];
  ifstream ifs(FileName);

  if (ifs.fail())
  {
    cerr << "Failed to open file." << endl;
    return -1;
  }

  int numInstanceFiles;
  int numTrial;
  string instance_file;
  int k;
  int mi;
  vector<string> InstanceFiles;
  vector<int> Ks;
  vector<int> maxIters;

  vector<vector<int> > Results;

  // ファイル読み込み
  ifs >> numInstanceFiles;      // 実行するインスタンスの数
  ifs >> numTrial;              // 1個のインスタンスを何回実行するか

  for (int i = 0; i < numInstanceFiles; i++)
  {
    ifs >> instance_file >> k >> mi;
    InstanceFiles.push_back(instance_file); // インスタンスのファイル
    Ks.push_back(k);                        // Kの値
    maxIters.push_back(mi);                 // 繰り返しの回数
  }


  // SCPのインスタンスを読み込む
  // char *FileName = argv[1];
  // FILE *SourceFile = fopen(FileName,"r");

  for (int i = 0; i < numInstanceFiles; i++)
  {
    string instance_file = InstanceFiles[i];
    int K = Ks[i];
    int maxIteration = maxIters[i];

    SCPinstance instance(instance_file);

    vector<int> result;


    for (int trial = 0; trial < numTrial; trial++)
    {
      // initialize
      for (int i = 0; i < instance.numColumns; i++) {
        SKCC.push_back(1);
        SCORE.push_back(0);
        TIMES.push_back(0);
      }
      for (int i = 0; i < instance.numRows; i++) {
        COST.push_back(1);
      }

      Rand rnd;
      //int seed = 0;
      rnd.seed(trial);
      // End Initialize;

      SCPsolution CSbest(instance, K);

      CSbest = DLL_com(instance, K, maxIteration, rnd);

      if (check_solution(instance, CSbest)) {
        //CSbest.print_solution();
        result.push_back(CSbest.totalWeight);
      }
    } // End trial

    Results.push_back(result);
  }

  // 出力
  for (int i = 0; i < numInstanceFiles; i++)
  {
    string instance_file = InstanceFiles[i];
    int K = Ks[i];
    int maxIteration = maxIters[i];

    cout << instance_file << "," << K << "," << maxIteration << ",";

    int Best_totalWeight = numeric_limits<int>::max();
    int Sum_totalWeight = 0;

    for (int t = 0; t < numTrial; t++)
    {
      if (Best_totalWeight > Results[i][t]) Best_totalWeight = Results[i][t];
      Sum_totalWeight += Results[i][t];
      cout << Results[i][t] << ",";
    }
    cout << Best_totalWeight << ","
         << (double)Sum_totalWeight / numTrial
         << endl;
  }


  return 0;
}
