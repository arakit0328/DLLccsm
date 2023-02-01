#include "SCPv.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

//
//
//  Class SCPinstance
//
//

// コンストラクタ
SCPinstance::SCPinstance(std::string instance_file)
{
  FILE *SourceFile = fopen(instance_file.c_str(), "r");

  if (SourceFile == NULL) throw DataException();
  else
  {
    int R, C, cost;
    fscanf(SourceFile, "%d", &R);
    fscanf(SourceFile, "%d", &C);
    numRows = R;
    numColumns = C;

    int* nCov = new int [numColumns];
    int* idx = new int [numColumns];
    for (int j = 0; j < numColumns; j++) {
      nCov[j] = 0;
      idx[j] = 0;
    }

    for (int j = 0; j < numColumns; j++)
    {
      std::vector<int> nn;
      Neighborhood.push_back(nn);
    }

    // read costs
    for(int j = 0; j < C; j++)
    {
      fscanf(SourceFile, "%d", &cost);
      Weight.push_back(cost);
    }

    // ファイルから各行の情報を読む
    int CoverNo, CoverID;
    std::vector<int> cov;


    for(int i = 0;  i < numRows; i++)
    {
      if (fscanf(SourceFile, "%d", &CoverNo) == EOF)
        throw (DataException());

      cov.clear();
      for(int j = 0; j < CoverNo; j++)
      {
        if (fscanf(SourceFile,"%d", &CoverID) == EOF)
          throw (DataException());

        if (CoverID >= 1 && CoverID <= numColumns)
        {
          cov.push_back(CoverID - 1);
          nCov[CoverID - 1]++;
        }
        else
          throw (DataException());
      }

      RowCovers.push_back(cov);
    }
    // ファイルの読み込み終了

    // 列の情報を作成
    for (int j = 0; j < numColumns; ++j) {
      std::vector<int> c(nCov[j]);
      ColEntries.push_back(c);
    }

    for (int i = 0; i < numRows; i++)
    {
      for (int c : RowCovers[i]) {
        ColEntries[c][idx[c]] = i;
        idx[c]++;
      }
    }
    // 列の情報の作成終了
    delete [] nCov;
    delete [] idx;
  }


  // 近傍を作る
  for (int j = 0; j < numColumns; ++j)
  {
    for (int r : ColEntries[j])
    {
      for (int c : RowCovers[r])
      {
        if (c != j) Neighborhood[j].push_back(c);
      }
    }
  }

  for (int j = 0; j < numColumns; ++j)
  {
    std::sort(Neighborhood[j].begin(), Neighborhood[j].end());
    Neighborhood[j].erase(std::unique(Neighborhood[j].begin(), Neighborhood[j].end()), Neighborhood[j].end());
  }
  // 近傍の生成終了

  

  
  // 簡単な方法でデータの正しさを確認
  long Sum1 = 0, Sum2 = 0;
  for (int i = 0; i < numRows; i++) Sum1 += RowCovers[i].size();
  for (int j = 0; j < numColumns; j++) Sum2 += ColEntries[j].size();

  //Incorrect Source File!
  if (Sum1 != Sum2)  throw (DataException());

  // 密度の計算
  Density = (float)Sum1/(numColumns * numRows);
}

// End: コンストラクタ


// デストラクタ
//SCPinstance::~SCPinstance()
// End: デストラクタ

// End SCPinstance

//
//
// Class SCPsolution
//
//

// コンストラクタ
SCPsolution::SCPsolution(const SCPinstance &inst, int k)
{
  //instance = inst;

  nRow = inst.numRows;
  nCol = inst.numColumns;
  K = k;
  num_Cover = 0;
  totalWeight = 0;

  for (int j = 0; j < nCol; ++j)
  {
    SOLUTION.push_back(0);
  }

  for (int i = 0; i < nRow; i++)
  {
    COVERED.push_back(0);
  }
}


// デストラクタ
SCPsolution::~SCPsolution()
{
}


// 候補解を初期化
void SCPsolution::initialize(SCPinstance &inst)
{
  num_Cover = 0;

  for (int j = 0; j < nCol; ++j)
  {
    SOLUTION[j] = 0;
  }

  for (int i = 0; i < nRow; ++i)
  {
    COVERED[i] = 0;
  }

  CS.clear();
}


// CSに列cを追加する
void SCPsolution::add_column(SCPinstance &inst, int c)
{
  if (SOLUTION[c])
  {
    printf("Column %d has already contained in CS\n", c);
    exit(1);
  }

  totalWeight += inst.Weight[c];

  SOLUTION[c] = 1;

  // CSに列cを追加
  CS.push_back(c);

  for (int r : inst.ColEntries[c])
  {
    COVERED[r]++;
    if (COVERED[r] == K)
    {
      num_Cover++;		// カバーされる行の数が増える
    }
  }
} // End add_column


// CSから列cを削除する
void SCPsolution::remove_column(SCPinstance &inst, int c)
{
  if (SOLUTION[c] == 0)
  {
    printf("Column %d is not contained in CS\n", c);
    exit(1);
  }

  totalWeight -= inst.Weight[c];

  SOLUTION[c] = 0;

  // CSから列cを削除
  CS.erase(remove(CS.begin(), CS.end(), c), CS.end());

  for (int r : inst.ColEntries[c])
  {
    COVERED[r]--;
    if (COVERED[r] == K-1)
    {
      num_Cover--;        // カバーされる行の数が減る
    }
  }
} // End remove_column


// CSの中身を表示
void SCPsolution::print_solution()
{
  sort(CS.begin(), CS.end());
  for (int c : CS)
  {
    printf("%d ", c + 1);
  }
  printf("\n");
} // End print_solution
