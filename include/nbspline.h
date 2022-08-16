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
            double rts; // time 
            double pos; // Position vector 1d
            double vel; // Velocity vector 1d
            double acc; // Acceleration vector 1d
        };

        struct nbs_pva_state_3d
        {
            vector<double> rts; // Relative time vector
            vector<Eigen::Vector3d> pos; // Position vector 3d
            vector<Eigen::Vector3d> vel; // Velocity vector 3d
            vector<Eigen::Vector3d> acc; // Acceleration vector 3d
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
        * @param t vector of time points, ti-order to ti the number of time points is k+1
        * etc 5 time points, 4 control points, 3rd order
        **/
        inline Eigen::MatrixXd create_general_m(int order, vector<double> t)
        {
            int k = order + 1;
            // According to the paper only able to calculate to order 3
            if (k > 4)
                k = 4;
            
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
                    // M(0,1) = 1 - M(0,0) - M(0,2);
                    M(0,1) = 4.0/6.0;
                    M(0,3) = 0.0;
                    
                    M(1,0) = -3 * M(0,0);
                    M(1,2) = 3 * ((t[3] - t[2]) * (t[2] - t[1])) / ((t[4] - t[1]) * (t[3] - t[1]));
                    M(1,1) = 3 * M(0,0) - M(1,2);
                    M(1,3) = 0.0;

                    M(2,0) = 3 * M(0,0);
                    M(2,2) = 3 * (t[3] - t[2]) / ((t[4] - t[1]) * (t[3] - t[1]));
                    M(2,1) = -3 * M(0,0) - M(2,2);
                    M(2,3) = 0.0;

                    M(3,3) = pow((t[3] - t[2]),2) / ((t[5] - t[2]) * (t[4] - t[2]));
                    M(3,2) = - M(2,2) / 3 - M(3,3) - pow((t[3] - t[2]),2) / ((t[4] - t[2]) * (t[4] - t[1]));
                    M(3,0) = - M(0,0);
                    M(3,1) = M(0,0) - M(3,2) - M(3,3);
                    break;
                default:
                    break;
            }

            return M;
        }

        bool check_query_time(vector<double> time, double query, std::pair<double,double>& t_i, int& off)
        {
            for (int i = 0; i < (int)time.size()-1; i++)
            {
                // Only considering [0,1) hence not including 1
                if (time[i+1] - query > 0 && query - time[i] >= 0)
                // if (time[i+1] - query > 0)
                {
                    t_i.first = time[i];
                    t_i.second = time[i+1];
                    off = i;
                    return true;
                }
            }

            return false;
        }

        /** @brief Create the pva state of a 1d non-uniform bspline 
         * @param time size is order + 2
         * @param cp size is order + 1
         * @param query_time double time that is queried
        **/
        inline nbs_pva_state_1d get_nbspline_1d(
            int order, vector<double> time, 
            vector<double> cp, double query_time)
        {
            nbs_pva_state_1d s;
            
            std::pair<double,double> t_i;
            int time_index_offset = 0;
            if (!check_query_time(time, query_time, t_i, time_index_offset))
                return s;
            std::cout << "time_index_offset {" << time_index_offset
                << "} t_i {" << t_i.first << " " << t_i.second << "}" << std::endl;

            int k = order + 1;
            vector<double> time_trim;
            for (int i = 0; i < k+1; i++)
                time_trim.push_back(time[time_index_offset+i]);

            Eigen::MatrixXd M = create_general_m(order, time_trim);
            
            // int n = (int)cp.size() - 1;
            // int m = n + order + 1; 

            // Only considering [0,1) hence not including 1
            double u_t = (query_time - t_i.first) / (t_i.second - t_i.first);

            // Control Points in a Span Column vector
            Eigen::VectorXd p = Eigen::VectorXd::Zero(k);
            // Position Row Vector, Velocity Row Vector, Acceleration Row Vector, Snap Row Vector
            Eigen::RowVectorXd u, du, ddu;
            u = du = ddu = Eigen::RowVectorXd::Zero(k); 

            // Make the u, du, ddu and p matrix
            for (int l = 0; l < k; l++)
            {
                u(l) = pow(u_t, l);
                p(l) = cp[time_index_offset + l];
                if (l >= 1)
                    du(l) = (l) * pow(u_t, l-1);
                if (l >= 2)
                    ddu(l) = (l) * (l-1) * pow(u_t, l-2);
            }

            s.pos = position_at_time_segment(M, u, p);
            s.vel = velocity_at_time_segment(M, du, p);
            s.acc = acceleration_at_time_segment(M, ddu, p);
            s.rts = query_time;

            return s;

        }

        /** @brief 
         * Calculate position value 
         * @param
         * M : order+1 x order+1 matrix
         * u : Row vector of position association vector
         * p : Control points in that segment
         * u * M * p  
        **/
        inline double position_at_time_segment(
            Eigen::MatrixXd M, Eigen::RowVectorXd u, Eigen::VectorXd p)
        {
            return (u * M * p)(0,0);
        }

        /** @brief 
         * Calculate velocity value
         * @param
         * M : order+1 x order+1 matrix
         * du : Row vector of velocity association vector
         * p : Control points in that segment
         * du * M * p  
        **/
        inline double velocity_at_time_segment(
            Eigen::MatrixXd M, Eigen::RowVectorXd du, Eigen::VectorXd p)
        {
            return (du * M * p)(0,0);
        }

        /** @brief 
         * Calculate acceleration value
         * @param
         * M : order+1 x order+1 matrix
         * ddu : Row vector of acceleration association vector
         * p : Control points in that segment
         * ddu * M * p  
        **/
        inline double acceleration_at_time_segment(
            Eigen::MatrixXd M, Eigen::RowVectorXd ddu, Eigen::VectorXd p)
        {
            return (ddu * M * p)(0,0);
        }

    };
}

#endif