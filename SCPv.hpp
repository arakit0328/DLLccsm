//---------------------------------------------------------------------------
// SCPインスタンスを扱うクラスの例
// vector
// 配列と構造体を使っています。
// クラスと言いながら，実際は構造体みたいな使い方
// Araki
//---------------------------------------------------------------------------
#pragma once

#include "Random.hpp"
#include <vector>
#include <string>
#include <cstdio>

//
//
//  Class SCPinstance  SCPのインスタンスを管理するクラス
//
//
class SCPinstance
{
public:
  int     numRows;             // The num of rows
  int     numColumns;          // The num of columns
  double  Density;             // The density of the matrix

public:
  SCPinstance(std::string SourceFile);
  ~SCPinstance() {}

  std::vector<std::vector<int> > RowCovers;	// 行をカバーする列のリスト
  std::vector<std::vector<int> > ColEntries;	// 列がカバーする行のリスト
  std::vector<std::vector<int> > Neighborhood;  // 列がカバーする行のリスト
  std::vector<int> Weight;                      // 列のコスト
};

class DataException {};


//
//
//  Class SCPsolution SCPの候補解を管理するクラス
//
//
class SCPsolution
{
public:
  //const SCPinstance& instance;        // インスタンス
  int nRow;                     // 行数
  int nCol;                     // 列数
  int K;                        // カバーする数
  int     totalWeight;

  std::vector<int> CS;                   // CS: 候補解（列番号のリスト）
  std::vector<int> SOLUTION;             // SOLUTION[j] = 1: 列jが候補解に含まれる
  std::vector<int> COVERED;              // COVERED[i]: 行iがカバーされている回数
  int num_Cover;                         // カバーされた行の数

public:
  // インスタンス，K, 行重みしきい値, oblivious_ratio
  SCPsolution(const SCPinstance &pData, int k);
  ~SCPsolution();

  // 候補解を初期化
  void initialize(SCPinstance &pData);

  // CSに含まれない列から最大スコアのものを選んで返す
  //int get_column_maxscore(SCPinstance &pData, Rand& rnd);

  // CSに含まれない列から最大スコアのものを選んで返す
  //int get_column_BMS(SCPinstance &pData, Rand& rnd);

  // CSに含まれない列から alpha * 最大スコア 以上のものを選んで返す
  //int get_column_grasp(SCPinstance &pData, double alpha, Rand& rnd);

  // CSに列cを追加する
  void add_column(SCPinstance &pData, int c);

  // CSから列cを削除する
  void remove_column(SCPinstance &pData, int c);

  // CSの中身を表示
  void print_solution();
};
