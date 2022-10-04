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
            int order, vector<t_p_sc> time, t_p_sc query, std::pair<t_p_sc, t_p_sc>& t_i, int& off)
        {
            for (int i = 0; i < (int)time.size()-1; i++)
            {
                double ms0 = duration<double>(query - time[i]).count();
                double ms1 = duration<double>(time[i+1] - query).count();
                // Only considering [0,1) hence not including 1
                if (ms0 >= 0 && ms1 > 0)
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

        /** @brief Create the pva state of a 1d non-uniform bspline from query time
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * @param order is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points 
         * @param query_time is the current point in time that is being queried
         * @param start is the start time 
        **/
        inline nbs_pva_state_1d get_nbspline_1d(
            int order, vector<t_p_sc> time, vector<double> cp, t_p_sc query_time, t_p_sc start)
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
            int k = order + 1;
            Eigen::MatrixXd M;
            double u_t, dt;
            int time_index_offset;

            if (!assemble_M_ut_dt_matrix(order, time, query_time, start, M, u_t, dt, time_index_offset))
                return s;

            // Control Points in a Span Column vector
            Eigen::VectorXd p = Eigen::VectorXd::Zero(k);
            // Position Row Vector, Velocity Row Vector, Acceleration Row Vector, Snap Row Vector
            Eigen::RowVectorXd u, du, ddu;
            assemble_u_p_matrix(u_t, time_index_offset, k, cp, u, du, ddu, p);

            s.pos = position_at_time_segment(u, M, p);
            s.vel = velocity_at_time_segment(dt, du, M, p);
            s.acc = acceleration_at_time_segment(dt, ddu, M, p);
            s.rts = query_time;

            return s;

        }

        /** @brief Create the pva state of a 1d non-uniform bspline from interval
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * @param order is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points 
         * @param interval is the interval used to move forward in time to find all points
         * @param start is the start time
        **/
        inline vector<nbs_pva_state_1d> get_nbspline_1d_all(
            int order, vector<t_p_sc> time, vector<double> cp, double interval, t_p_sc start)
        {
            vector<nbs_pva_state_1d> v_s;
            // According to the paper only able to calculate to order 3
            if (order > 3)
                return v_s;
            if (cp.empty() || time.empty())
                return v_s;
            if (time.size() != cp.size() + (order-1))
            {
                std::cout << KRED << "time vector size not correct!" << KNRM << std::endl;
                return v_s;
            }

            double total_duration = duration<double>(time.back() - time.front()).count();
            int divisions = (int)floor(total_duration / interval);

            for (int i = 1; i < divisions; i++)
            {
                int k = order + 1;
                Eigen::MatrixXd M;
                double u_t, dt;
                int time_index_offset;

                t_p_sc query_time = time.front() + milliseconds(i * (int)round(interval*1000));
                
                nbs_pva_state_1d s;

                if (!assemble_M_ut_dt_matrix(order, time, query_time, start, M, u_t, dt, time_index_offset))
                    return v_s;

                // Control Points in a Span Column vector
                Eigen::VectorXd p = Eigen::VectorXd::Zero(k);
                // Position Row Vector, Velocity Row Vector, Acceleration Row Vector, Snap Row Vector
                Eigen::RowVectorXd u, du, ddu;
                assemble_u_p_matrix(u_t, time_index_offset, k, cp, u, du, ddu, p);

                s.pos = position_at_time_segment(u, M, p);
                s.vel = velocity_at_time_segment(dt, du, M, p);
                s.acc = acceleration_at_time_segment(dt, ddu, M, p);
                s.rts = query_time;
            
                v_s.push_back(s);
            }

            return v_s;
        }

        /** @brief Create the pva state of a 1d non-uniform bspline
         * This is with prior calculation of M matrix and u_t and offset
         * @param dt is the interval for the current pair of knots
         * @param time_index_offset is the time_index_offset which is the (floor)time_point index that current time is in
         * @param cp is the vector of control points 
         * @param u_t is current knot factor
         * @param M is basis matrix
         * @param k is order + 1
        **/
        inline nbs_pva_state_1d get_nbspline_1d_w_prior(
            double dt, int time_index_offset, vector<double> cp, double u_t, Eigen::MatrixXd M, int k)
        {
            nbs_pva_state_1d s;

            // Control Points in a Span Column vector
            Eigen::VectorXd p = Eigen::VectorXd::Zero(k);
            // Position Row Vector, Velocity Row Vector, Acceleration Row Vector, Snap Row Vector
            Eigen::RowVectorXd u, du, ddu;
            assemble_u_p_matrix(u_t, time_index_offset, k, cp, u, du, ddu, p);

            s.pos = position_at_time_segment(u, M, p);
            s.vel = velocity_at_time_segment(dt, du, M, p);
            s.acc = acceleration_at_time_segment(dt, ddu, M, p);

            return s;
        }

        /** @brief Create the pva state of a 3d non-uniform bspline from a query time
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * Using get_nbspline_1d_w_prior() to save repetition in computation of prior M matrix, u_t and offset 
         * @param order is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points in the 3d coordinates
         * @param query_time is the current point in time that is being queried 
         * @param start is the start time
        **/
        inline nbs_pva_state_3d get_nbspline_3d(
            int order, vector<t_p_sc> time, vector<Eigen::Vector3d> cp, t_p_sc query_time, t_p_sc start)
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
            int k = order + 1;
            Eigen::MatrixXd M;
            double u_t, dt;
            int time_index_offset;

            if (!assemble_M_ut_dt_matrix(order, time, query_time, start, M, u_t, dt, time_index_offset))
                return ss;

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

        /** @brief Create the pva state of a 3d non-uniform bspline from an interval
         * This is 1 pass function, does not calculate more that 1 instance in the bspline 
         * Using get_nbspline_1d_w_prior() to save repetition in computation of prior M matrix, u_t and offset 
         * @param order is degree of the spline
         * @param time is the knot vector that was acquired initially
         * @param cp is the vector of control points in the 3d coordinates
         * @param interval is the interval used to move forward in time to find all points
         * @param start is the start time
        **/
        inline vector<nbs_pva_state_3d> get_nbspline_3d_all(
            int order, vector<t_p_sc> time, vector<Eigen::Vector3d> cp, double interval, t_p_sc start)
        {
            vector<nbs_pva_state_3d> v_ss;

            // According to the paper only able to calculate to order 3
            if (order > 3)
                return v_ss;
            if (cp.empty() || time.empty())
                return v_ss;
            if (time.size() != cp.size() + (order-1))
            {
                std::cout << KRED << "time vector size not correct!" << KNRM << std::endl;
                return v_ss;
            }
            int k = order + 1;

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

            double total_duration = duration<double>(time.back() - time.front()).count();
            int divisions = (int)floor(total_duration / interval);

            for (int i = 1; i < divisions; i++)
            {
                nbs_pva_state_3d ss;
                Eigen::MatrixXd M;
                double u_t, dt;
                int time_index_offset;

                t_p_sc query_time = time.front() + milliseconds(i * (int)round(interval*1000));

                if (!assemble_M_ut_dt_matrix(order, time, query_time, start, M, u_t, dt, time_index_offset))
                    return v_ss;
                
                ss.rts = query_time;
                nbs_pva_state_1d x = get_nbspline_1d_w_prior(
                    dt, time_index_offset, rv.xcp, u_t, M, k);
                nbs_pva_state_1d y = get_nbspline_1d_w_prior(
                    dt, time_index_offset, rv.ycp, u_t, M, k);
                nbs_pva_state_1d z = get_nbspline_1d_w_prior(
                    dt, time_index_offset, rv.zcp, u_t, M, k);

                ss.pos = Eigen::Vector3d(x.pos, y.pos, z.pos);
                ss.vel = Eigen::Vector3d(x.vel, y.vel, z.vel);
                ss.acc = Eigen::Vector3d(x.acc, y.acc, z.acc);

                v_ss.push_back(ss);
            }

            return v_ss;
        }

        /** @brief 
         * Assemble the knot (u) and cp (p) vector for matrix multiplication 
         * @param u_t is current knot factor
         * @param t_i_o is the time_index_offset which is the (floor)time_point index that current time is in
         * @param k is order + 1 
         * @param u (return) the row vector of position association vector
         * @param p (return) the control points in that segment
        **/
        inline void assemble_u_p_matrix(
            double u_t, int t_i_o, int k, vector<double> cp, 
            Eigen::RowVectorXd& u, Eigen::RowVectorXd& du, Eigen::RowVectorXd& ddu, Eigen::VectorXd& p)
        {
            // Zero the values of u, du, ddu
            u = du = ddu = Eigen::RowVectorXd::Zero(k); 

            // Make the u, du, ddu and p matrix
            // std::cout << "P vector";
            for (int l = 0; l < k; l++)
            {
                u(l) = pow(u_t, l);
                p(l) = cp[t_i_o + l];
                // std::cout << " " << p(l);
                if (l >= 1)
                    du(l) = (l) * pow(u_t, l-1);
                if (l >= 2)
                    ddu(l) = (l) * (l-1) * pow(u_t, l-2);
            }
            // std::cout << std::endl;
        }

        /** @brief 
         * Assemble the M matrix, u_t (factor), and dt
         * @param o is the degree/order
         * @param t is the knot vector
         * @param q is the current time
         * @param s is the start time
         * @param M (return) the general matrix for the basis
         * @param u_t (return) current knot factor
         * @param dt (return) the difference between the knots that the current point is within
         * @param t_i_o (return) the time_index_offset which is the (floor)time_point index that current time is in
        **/
        inline bool assemble_M_ut_dt_matrix(
            int o, vector<t_p_sc> t, t_p_sc q, t_p_sc s, 
            Eigen::MatrixXd& M, double& u_t, double& dt, int& t_i_o)
        {
            int k = o + 1;
            std::pair<t_p_sc, t_p_sc> t_i;
            if (!check_query_time(o, t, q, t_i, t_i_o))
                return false;

            vector<double> t_t;
            // std::cout << "time_trim vector";
            for (int i = 0; i < k + (o - 1); i++)
            {
                double specific_time = duration<double>(t[t_i_o + i] - s).count();
                t_t.push_back(specific_time);
                // std::cout << " " << specific_time;
            }
            // std::cout << std::endl;

            M = create_general_m(o, t_t);

            // u_t = (query_time - t_i.first) / (t_i.second - t_i.first)
            double numerator = duration<double>(q - t_i.first).count();
            double denominator = duration<double>(t_i.second - t_i.first).count();
            // Only considering [0,1) hence not including 1
            u_t = numerator / denominator;
            dt = denominator;

            return true;
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