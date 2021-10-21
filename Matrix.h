#pragma once

#include <iostream>

namespace custom
{
/**
 * @brief 矩阵
 *
 */
class  matrix
{
public:
    matrix(unsigned int m = 0, unsigned int n = 0);
    matrix(const matrix& mat);
    ~matrix();
    /**
     * @brief 重新调整大小
     *
     * @param m
     * @param n
     */
    void resize(unsigned int m, unsigned int n);
    /**
     * @brief 清除
     *
     */
    void clear();
    /**
     * @brief 设为单位阵
     *
     */
    void load_identity();
    /**
    * 小心使用，参数必须显式指定为double型，不可使用任何自动类型转换
    */
    void init(double d, ...);

    /**
     * @brief 列数
     *
     * @return unsigned int
     */
    unsigned int get_columns(void) const;
    /**
     * @brief 行数
     *
     * @return unsigned int
     */
    unsigned int get_lines(void) const;
    /**
     * @brief 取得数据
     *
     * @return const double
     */
    const double* get_data(void) const;

    /**
     * @brief 取得元素
     *
     * @param i 行号，从0开始
     * @param j 列号，从0开始
     * @return double
     */
    double at(unsigned i, unsigned j) const;
    /**
     * @brief 取得元素
     * @see at(unsigned,unsigned) const
     * @param i
     * @param j
     * @return double 元素引用
     */
    double& at(unsigned int i, unsigned int j);
    /**
     * @brief 取得元素
     * @see at(unsigned int, unsigned int)
     * @param i
     * @param j
     * @return double &operator
     */
    double& operator()(unsigned int i, unsigned int j);
    /**
     * @brief 取得行
     *
     * @param i
     * @return double *operator
     */
    double* operator[](unsigned int i);
    matrix operator*(const matrix& A) const;
    matrix operator*(double a);
    matrix operator+(matrix& mat);
    matrix& operator+=(matrix& mat);
    matrix operator-(matrix& mat);
    matrix operator-();
    matrix& operator-=(matrix& mat);
    matrix& operator=(const matrix& A);
    matrix& operator=(double* data/* line first */);

    friend std::ostream& operator<<(std::ostream& out, matrix& mat);

    /**
     * @brief 求逆
     *
     * @return matrix
     */
    matrix inverse();
    /**
     * @brief inverse4X4
     * 适用于 4X4齐次正交矩阵 求逆
     * @return
     */
    matrix inverse4X4();
    /**
     * @brief 转置
     *
     * @return matrix
     */
    matrix transpose();
    /**
     * @brief Hermite标准型
     *
     * @return matrix
     */
    matrix Hermite_canonical_form();
    /**
     * @brief 伪逆
     *
     * @return matrix
     */
    matrix pseudo_inverse();
    /**
     * @brief 为空
     *
     * @return bool
     */
    bool null();
    /**
     * @brief 数据精度
     *
     * @param pc
     */
    void set_precision(double pc);
    /**
     * @brief 满秩分解
     *
     * @param B
     * @param L
     * @param R
     */
    static void FullRankDecomposition( matrix& B, matrix& L, matrix& R );

    /**
     * @brief setrowVal 设置每行数据
     * @param row
     * @param vrow
     */
    void setrowVal(int row, float vrow[]);

    /**
     * @brief setValue 设置单个元素数值
     * @param row 行
     * @param col 列
     * @param val 值
     */
    void  setValue(int row, int col, float val);
    /**
     * @brief getValue 得到元素值
     * @param row 行
     * @param col 列
     * @return 值
     */
    float getValue(int row, int col);

protected:
    double* data;
    unsigned int columns,lines;
    double eps/* = 1e-10 */;
    unsigned int rank;
};
}
