/*
* nbspline.hpp
*
* ---------------------------------------------------------------------
* Copyright (C) 2022 Matthew (matthewoots at gmail.com)
*
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License
*  as published by the Free Software Foundation; either version 2
*  of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
* ---------------------------------------------------------------------
*/

#ifndef NBSPLINE_H
#define NBSPLINE_H

#include <iostream>
#include <string>
#include <chrono>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

using namespace Eigen;
using namespace std;
using namespace std::chrono;

namespace nbspline
{
    typedef time_point<std::chrono::system_clock> t_p_sc; // giving a typename

    class bspline_trajectory
    {
        private:
        
        /*
        * General matrix representations for B-splines See Theorem 1 of page 182 for getting M
        * https://link.springer.com/article/10.1007/s003710050206
        */

        /** @brief 
         * Calculate position value 
         * @param u is the row vector of position association vector
         * @param M is the degree+1 x degree+1 matrix
         * @param p is the control points in that segment
         * u * M * p  
        **/
        inline double position_at_time_segment(
            Eigen::RowVectorXd u, Eigen::MatrixXd M, Eigen::VectorXd p)
        {
            return (u * M * p)(0,0);
        }

        /** @brief 
         * Calculate velocity value
         * @param dt knot interval
         * @param du is the row vector of velocity association vector
         * @param M is the degree+1 x degree+1 matrix
         * @param p is the control points in that segment
         * (1/dt) * du * M * p  
        **/
        inline double velocity_at_time_segment(
            double dt, Eigen::RowVectorXd du, Eigen::MatrixXd M, Eigen::VectorXd p)
        {
            return (1/dt) * (du * M * p)(0,0);
        }

        /** @brief 
         * Calculate acceleration value
         * @param dt knot interval
         * @param ddu is the row vector of acceleration association vector
         * @param M is the degree+1 x degree+1 matrix
         * @param p is the control points in that segment
         * pow((1/dt),2) * ddu * M * p  
        **/
        inline double acceleration_at_time_segment(
            double dt, Eigen::RowVectorXd ddu, Eigen::MatrixXd M, Eigen::VectorXd p)
        {
            return pow((1/dt),2) * (ddu * M * p)(0,0);
        }

        /** @brief Create the range of values within the array of index **/
        inline vector<int> int_range_to_vector(int min, int max)
        {
            vector<int> v;
            for (int i = min; i <= max; i++)
                v.push_back(i);
            
            return v;
        }

        public:

        struct nbs_pva_state_1d
        {
            t_p_sc rts; // time 
            double pos; // position vector 1d
            double vel; // velocity vector 1d
            double acc; // acceleration vector 1d
        };

        struct nbs_pva_state_3d
        {
            t_p_sc rts; // time
            Eigen::Vector3d pos; // position 3d
            Eigen::Vector3d vel; // velocity 3d
            Eigen::Vector3d acc; // acceleration 3d
        };

        struct row_vector_3d
        {
            vector<double> xcp; // x control point
            vector<double> ycp; // y control point
            vector<double> zcp; // z control point
        };

        /** @brief Creating Non-uniform Bspline basis M matrix
        * https://xiaoxingchen.github.io/2020/03/02/bspline_in_so3/general_matrix_representation_for_bsplines.pdf
        * @param degree the degree of the spline 
        * @param t vector of time points, t_{i-degree} to t_i the number of time points is k+1
        * etc 5 time points, 4 control points, 3rd degree
        **/
        Eigen::MatrixXd create_general_m(int degree, vector<double> t);

        /** @brief Find the current pair of knots that we are inbetween of
         * @param degree is used to lower the offset so that we can construct the M matrix according to Kaihuai Qin's formulation
         * @param time is the knot vector that was acquired initially
         * @param query is the current point in time that is being queried 
         * @param t_i (return) the current pair of knots
         * @param off (return) the current offset
        **/
        bool check_query_time(
            int degree, vector<t_p_sc> time, t_p_sc query, std::pair<t_p_sc, t_p_sc>& t_i, int& off);

        /** @brief Create the pva state of a 1d non-uniform bspline from query time
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * @param degree is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points 
         * @param query_time is the current point in time that is being queried
         * @param start is the start time 
        **/
        nbs_pva_state_1d get_nbspline_1d(
            int degree, vector<t_p_sc> time, vector<double> cp, t_p_sc query_time, t_p_sc start);

        /** @brief Create the pva state of a 1d non-uniform bspline from interval
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * @param order is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points 
         * @param interval is the interval used to move forward in time to find all points
         * @param start is the start time
        **/
        vector<nbs_pva_state_1d> get_nbspline_1d_all(
            int degree, vector<t_p_sc> time, vector<double> cp, double interval, t_p_sc start);

        /** @brief Create the pva state of a 1d non-uniform bspline
         * This is with prior calculation of M matrix and u_t and offset
         * @param dt is the interval for the current pair of knots
         * @param time_index_offset is the time_index_offset which is the (floor)time_point index that current time is in
         * @param cp is the vector of control points 
         * @param u_t is current knot factor
         * @param M is basis matrix
         * @param k is degree + 1
        **/
        nbs_pva_state_1d get_nbspline_1d_w_prior(
            double dt, int time_index_offset, vector<double> cp, double u_t, Eigen::MatrixXd M, int k);

        /** @brief Create the pva state of a 3d non-uniform bspline from a query time
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * Using get_nbspline_1d_w_prior() to save repetition in computation of prior M matrix, u_t and offset 
         * @param degree is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points in the 3d coordinates
         * @param query_time is the current point in time that is being queried 
         * @param start is the start time
        **/
        nbs_pva_state_3d get_nbspline_3d(
            int degree, vector<t_p_sc> time, vector<Eigen::Vector3d> cp, t_p_sc query_time, t_p_sc start);

        /** @brief Create the pva state of a 3d non-uniform bspline from an interval
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * Using get_nbspline_1d_w_prior() to save repetition in computation of prior M matrix, u_t and offset 
         * @param degree is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points in the 3d coordinates
         * @param interval is the interval used to move forward in time to find all points
         * @param start is the start time
        **/
        vector<nbs_pva_state_3d> get_nbspline_3d_all(
            int degree, vector<t_p_sc> time, vector<Eigen::Vector3d> cp, double interval, t_p_sc start);

        /** @brief 
         * Assemble the knot (u) and cp (p) vector for matrix multiplication 
         * @param u_t is current knot factor
         * @param t_i_o is the time_index_offset which is the (floor)time_point index that current time is in
         * @param k is degree + 1 
         * @param u (return) the row vector of position association vector
         * @param p (return) the control points in that segment
        **/
        void assemble_u_p_matrix(
            double u_t, int t_i_o, int k, vector<double> cp, 
            Eigen::RowVectorXd& u, Eigen::RowVectorXd& du, Eigen::RowVectorXd& ddu, Eigen::VectorXd& p);

        /** @brief 
         * Assemble the M matrix, u_t (factor), and dt
         * @param d is the degree
         * @param t is the knot vector
         * @param q is the current time
         * @param s is the start time
         * @param M (return) the general matrix for the basis
         * @param u_t (return) current knot factor
         * @param dt (return) the difference between the knots that the current point is within
         * @param t_i_o (return) the time_index_offset which is the (floor)time_point index that current time is in
        **/
        bool assemble_M_ut_dt_matrix(
            int d, vector<t_p_sc> t, t_p_sc q, t_p_sc s, 
            Eigen::MatrixXd& M, double& u_t, double& dt, int& t_i_o);

    };
}

#endif