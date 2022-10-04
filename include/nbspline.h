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

        /** @brief Create the range of values within the array of index **/
        inline vector<int> int_range_to_vector(int min, int max)
        {
            vector<int> v;
            for (int i = min; i <= max; i++)
                v.push_back(i);
            
            return v;
        }

        /** @brief Creating Non-uniform Bspline basis M matrix
        * https://xiaoxingchen.github.io/2020/03/02/bspline_in_so3/general_matrix_representation_for_bsplines.pdf
        * @param order the degree of the spline 
        * @param t vector of time points, ti-order to ti the number of time points is k+1
        * etc 5 time points, 4 control points, 3rd order
        **/
        inline Eigen::MatrixXd create_general_m(int order, vector<double> t)
        {
            int k = order + 1;
            
            Eigen::MatrixXd M = Eigen::MatrixXd::Zero(k,k);

            // Eigen (row,column)
            switch(k)
            {
                case 2:
                    M(0,0) = 1;
                    M(0,1) = 0;

                    M(1,0) = -1;
                    M(1,1) = 1;
                    break;
                case 3:
                    M(0,0) = (t[2] - t[1]) / (t[2] - t[0]);
                    M(0,1) = (t[1] - t[0]) / (t[2] - t[0]);
                    M(0,2) = 0.0;
                    
                    M(1,0) = -2 * (t[2] - t[1]) / (t[2] - t[0]);
                    M(1,1) = 2 * (t[2] - t[1]) / (t[2] - t[0]);
                    M(1,2) = 0.0;
                    
                    M(2,0) = (t[2] - t[1]) / (t[2] - t[0]);
                    M(2,1) = - (t[2] - t[1]) * (1/(t[2] - t[0]) + 1/(t[3] - t[1]));
                    // M(2,2) = (t[2] - t[1]) / (t[3] - t[0]);
                    M(2,2) = (t[2] - t[1]) / (t[2] - t[0]);
                    break;
                case 4:
                    M(0,0) = pow((t[3] - t[2]),2) / ((t[3] - t[1]) * (t[3] - t[0]));
                    M(0,2) = pow((t[2] - t[1]),2) / ((t[4] - t[1]) * (t[3] - t[1]));
                    M(0,1) = 1 - M(0,0) - M(0,2);
                    M(0,3) = 0.0;
                    
                    M(1,0) = -3 * M(0,0);
                    M(1,2) = 3 * ((t[3] - t[2]) * (t[2] - t[1])) / ((t[4] - t[1]) * (t[3] - t[1]));
                    M(1,1) = 3 * M(0,0) - M(1,2);
                    M(1,3) = 0.0;

                    M(2,0) = 3 * M(0,0);
                    M(2,2) = 3 * pow((t[3] - t[2]),2) / ((t[4] - t[1]) * (t[3] - t[1]));
                    M(2,1) = -3 * M(0,0) - M(2,2);
                    M(2,3) = 0.0;

                    M(3,3) = pow((t[3] - t[2]),2) / ((t[5] - t[2]) * (t[4] - t[2]));
                    M(3,2) = -(M(2,2) / 3) - M(3,3) - pow((t[3] - t[2]),2) / ((t[4] - t[2]) * (t[4] - t[1]));
                    M(3,0) = - M(0,0);
                    M(3,1) = M(0,0) - M(3,2) - M(3,3);
                    break;
                default:
                    break;
            }

            return M;
        }

        /** @brief Find the current pair of knots that we are inbetween of
         * @param order is used to lower the offset so that we can construct the M matrix according to Kaihuai Qin's formulation
         * @param time is the knot vector that was acquired initially
         * @param query is the current point in time that is being queried 
         * @param t_i (return) the current pair of knots
         * @param off (return) the current offset
        **/
        bool check_query_time(
            int order, vector<t_p_sc> time, t_p_sc query, 
            std::pair<t_p_sc, t_p_sc>& t_i, int& off)
        {
            for (int i = 0; i < (int)time.size()-1; i++)
            {
                auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(query - time[i]);
                auto ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(time[i+1] - query);
                // Only considering [0,1) hence not including 1
                // Change it from milliseconds to seconds
                if ((double)ms0.count()/1000.0 >= 0 && (double)ms1.count()/1000.0 > 0)
                {
                    t_i.first = time[i];
                    t_i.second = time[i+1];
                    // off = i;
                    off = i-(order-1);
                    return true;
                }
            }

            return false;
        }

        /** @brief Create the pva state of a 1d non-uniform bspline
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * @param order is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points 
         * @param query_time is the current point in time that is being queried 
        **/
        inline nbs_pva_state_1d get_nbspline_1d(
            int order, vector<t_p_sc> time, 
            vector<double> cp, t_p_sc query_time, t_p_sc start)
        {
            nbs_pva_state_1d s;
            // According to the paper only able to calculate to order 3
            if (order > 3)
                return s;
            if (cp.empty() || time.empty())
                return s;
            if (time.size() != cp.size() + (order-1))
            {
                std::cout << KRED << "time vector size not correct!" << KNRM << std::endl;
                return s;
            }
            
            std::pair<t_p_sc, t_p_sc> t_i;
            int time_index_offset = 0;
            if (!check_query_time(order, time, query_time, t_i, time_index_offset))
                return s;

            int k = order + 1;
            vector<double> time_trim;
            // std::cout << "time_trim vector";
            for (int i = 0; i < k+(order-1); i++)
            {
                double specific_time = duration<double>(time[time_index_offset+i] - start).count();
                time_trim.push_back(specific_time);
                // std::cout << " " << specific_time;
            }
            // std::cout << std::endl;

            Eigen::MatrixXd M = create_general_m(order, time_trim);

            // u_t = (query_time - t_i.first) / (t_i.second - t_i.first)
            double numerator = duration<double>(query_time - t_i.first).count();
            double denominator = duration<double>(t_i.second - t_i.first).count();
            // Only considering [0,1) hence not including 1
            double u_t = numerator / denominator;
            double dt = denominator;

            // Control Points in a Span Column vector
            Eigen::VectorXd p = Eigen::VectorXd::Zero(k);
            // Position Row Vector, Velocity Row Vector, Acceleration Row Vector, Snap Row Vector
            Eigen::RowVectorXd u, du, ddu;
            u = du = ddu = Eigen::RowVectorXd::Zero(k); 

            // Make the u, du, ddu and p matrix
            // std::cout << "P vector";
            for (int l = 0; l < k; l++)
            {
                u(l) = pow(u_t, l);
                p(l) = cp[time_index_offset + l];
                // std::cout << " " << p(l);
                if (l >= 1)
                    du(l) = (l) * pow(u_t, l-1);
                if (l >= 2)
                    ddu(l) = (l) * (l-1) * pow(u_t, l-2);
            }
            // std::cout << std::endl;

            s.pos = position_at_time_segment(u, M, p);
            s.vel = velocity_at_time_segment(dt, du, M, p);
            s.acc = acceleration_at_time_segment(dt, ddu, M, p);
            s.rts = query_time;

            return s;

        }

        /** @brief Create the pva state of a 1d non-uniform bspline
         * This is with prior calculation of M matrix and u_t and offset
         * @param order is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points 
         * @param query_time is the current point in time that is being queried 
        **/
        inline nbs_pva_state_1d get_nbspline_1d_w_prior(
            double dt, int time_index_offset, vector<double> cp, 
            double u_t, Eigen::MatrixXd M, int k)
        {
            nbs_pva_state_1d s;

            // Control Points in a Span Column vector
            Eigen::VectorXd p = Eigen::VectorXd::Zero(k);
            // Position Row Vector, Velocity Row Vector, Acceleration Row Vector, Snap Row Vector
            Eigen::RowVectorXd u, du, ddu;
            u = du = ddu = Eigen::RowVectorXd::Zero(k); 

            // Make the u, du, ddu and p matrix
            // std::cout << "P vector";
            for (int l = 0; l < k; l++)
            {
                u(l) = pow(u_t, l);
                p(l) = cp[time_index_offset + l];
                // std::cout << " " << p(l);
                if (l >= 1)
                    du(l) = (l) * pow(u_t, l-1);
                if (l >= 2)
                    ddu(l) = (l) * (l-1) * pow(u_t, l-2);
            }
            // std::cout << std::endl;

            s.pos = position_at_time_segment(u, M, p);
            s.vel = velocity_at_time_segment(dt, du, M, p);
            s.acc = acceleration_at_time_segment(dt, ddu, M, p);

            return s;

        }

        /** @brief Create the pva state of a 3d non-uniform bspline
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * Using get_nbspline_1d_w_prior() to save repetition in computation of prior M matrix, u_t and offset 
         * @param order is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points in the 3d coordinates
         * @param query_time is the current point in time that is being queried 
        **/
        inline nbs_pva_state_3d get_nbspline_3d(
            int order, vector<t_p_sc> time, vector<Eigen::Vector3d> cp, 
            t_p_sc query_time, t_p_sc start)
        {
            nbs_pva_state_3d ss;
            ss.rts = query_time;

            // According to the paper only able to calculate to order 3
            if (order > 3)
                return ss;
            if (cp.empty() || time.empty())
                return ss;
            if (time.size() != cp.size() + (order-1))
            {
                std::cout << KRED << "time vector size not correct!" << KNRM << std::endl;
                return ss;
            }
            
            std::pair<t_p_sc, t_p_sc> t_i;
            int time_index_offset = 0;
            if (!check_query_time(order, time, query_time, t_i, time_index_offset))
                return ss;

            int k = order + 1;
            vector<double> time_trim;
            // std::cout << "time_trim vector";
            for (int i = 0; i < k+(order-1); i++)
            {
                double specific_time = duration<double>(time[time_index_offset+i] - start).count();
                time_trim.push_back(specific_time);
                // std::cout << " " << specific_time;
            }
            // std::cout << std::endl;

            Eigen::MatrixXd M = create_general_m(order, time_trim);

            // u_t = (query_time - t_i.first) / (t_i.second - t_i.first)
            double numerator = duration<double>(query_time - t_i.first).count();
            double denominator = duration<double>(t_i.second - t_i.first).count();
            // Only considering [0,1) hence not including 1
            double u_t = numerator / denominator;
            double dt = denominator;

            row_vector_3d rv;
            // time_point<std::chrono::system_clock> t_s = system_clock::now();
            // Reorganize the control points into vectors that can be passed into 1d_bspline function
            for (int i = 0; i < (int)cp.size(); i++)
            {
                rv.xcp.push_back(cp[i].x());
                rv.ycp.push_back(cp[i].y());
                rv.zcp.push_back(cp[i].z());
            }
            // auto t_c = duration<double>(system_clock::now() - t_s).count() * 1000;
            // std::cout << "conversion time " << KGRN << t_c << "ms" << KNRM << std::endl;


            nbs_pva_state_1d x = get_nbspline_1d_w_prior(
                dt, time_index_offset, rv.xcp, u_t, M, k);
            nbs_pva_state_1d y = get_nbspline_1d_w_prior(
                dt, time_index_offset, rv.ycp, u_t, M, k);
            nbs_pva_state_1d z = get_nbspline_1d_w_prior(
                dt, time_index_offset, rv.zcp, u_t, M, k);

            ss.pos = Eigen::Vector3d(x.pos, y.pos, z.pos);
            ss.vel = Eigen::Vector3d(x.vel, y.vel, z.vel);
            ss.acc = Eigen::Vector3d(x.acc, y.acc, z.acc);

            return ss;
        }

        /** @brief 
         * Calculate position value 
         * @param u is the row vector of position association vector
         * @param M is the order+1 x order+1 matrix
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
         * @param M is the order+1 x order+1 matrix
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
         * @param M is the order+1 x order+1 matrix
         * @param p is the control points in that segment
         * pow((1/dt),2) * ddu * M * p  
        **/
        inline double acceleration_at_time_segment(
            double dt, Eigen::RowVectorXd ddu, Eigen::MatrixXd M, Eigen::VectorXd p)
        {
            return pow((1/dt),2) * (ddu * M * p)(0,0);
        }

    };
}

#endif