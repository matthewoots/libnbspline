/*
* nbspline.cpp
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

#include "nbspline.h"

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
    bool bspline_trajectory::find_nearest_knot(
        t_p_sc tp, vector<t_p_sc> v, t_p_sc &found, int &index)
    {
        for (int i = 0; i < (int)v.size(); i++)
        {
            if (duration<double>(tp - v[i]).count() < 0)
            {
                if ((i-1) >= 0)
                {
                    found = v[i-1];
                    index = i-1;
                    return true;
                }
                else
                    return false;
            }
        }

        return false;
    }

    void bspline_trajectory::distribute_3d_control_points(
        int degree, Eigen::Vector3d p, vector<Eigen::Vector3d> wp, 
        double m_v, double k_s, vector<double> &time_waypoint,
        vector<Eigen::Vector3d> &control_points)
    {
        double est_dist_knot = m_v * k_s;
        // std::cout << "est_dist_knot = " << est_dist_knot << std::endl;

        // Add current position into the waypoint list
        vector<Eigen::Vector3d> tmp = wp;
        tmp.insert(tmp.begin(), p);
        
        vector<double> seg_dist;
        vector<int> seg_size;
        vector<Eigen::Vector3d> seg_vect;
        int total_segment = 0;
        
        for (int i = 0; i < (int)tmp.size() - 1; i++)
        {
            seg_dist.push_back((tmp[i+1] - tmp[i]).norm());
            seg_vect.push_back((tmp[i+1] - tmp[i]).normalized());
        }
        
        int number_of_segments = (int)seg_dist.size();

        for (int i = 0; i < number_of_segments; i++)
        {
            // ceil helps to push values above 0 to 1 or more
            // or else segment count is 0 and causes an error
            int knots_in_seg = (int)floor(seg_dist[i]/est_dist_knot);
            total_segment += knots_in_seg;
            
            seg_size.push_back(knots_in_seg);
        }

        control_points.push_back(p);
        time_waypoint.push_back(0.0);
        int count = 0;
        for (int i = 0; i < number_of_segments; i++)
        {
            if (seg_size[i] > 0)
            {
                for (int j = 0; j < seg_size[i]; j++)
                {
                    count += 1;
                    time_waypoint.push_back(count * k_s);
                    control_points.push_back(
                        tmp[i] + (j+1) * seg_vect[i] * est_dist_knot);
                }
            }
            
            count += 1;
            time_waypoint.push_back(count * k_s);
            control_points.push_back(tmp[i+1]);
        }

        // for (int i = 0; i < (int)control_points.size(); i++)
        //     std::cout << control_points[i].transpose() << " ";
        // std:cout << std::endl;

        // Since the relationship of time_point_size and cp_size is cp_size + (degree-1)
        // for (int i = 0; i < (degree-1); i++)
        // {
        //     int last_index = time_waypoint.size()-1;
        //     time_waypoint.push_back(last_index + k_s);
        // }

        return;
    }

    Eigen::MatrixXd bspline_trajectory::create_general_m(int degree, vector<double> t)
    {
        int k = degree + 1;
        
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

    bool bspline_trajectory::check_query_time(
        int degree, vector<t_p_sc> time, t_p_sc query, std::pair<t_p_sc, t_p_sc>& t_i, int& off)
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
                off = i-(degree-1);
                return true;
            }
        }

        return false;
    }

    bspline_trajectory::nbs_pva_state_1d bspline_trajectory::get_nbspline_1d(
        int degree, vector<t_p_sc> time, vector<double> cp, t_p_sc query_time, t_p_sc start)
    {
        bspline_trajectory::nbs_pva_state_1d s;
        // According to the paper only able to calculate to degree 3 which is order 4
        if (degree > 3)
            return s;
        if (cp.empty() || time.empty())
            return s;
        if (time.size() != cp.size() + (degree-1))
        {
            std::cout << KRED << "time vector size not correct!" << KNRM << std::endl;
            return s;
        }
        int k = degree + 1;
        Eigen::MatrixXd M;
        double u_t, dt;
        int time_index_offset;

        if (!assemble_M_ut_dt_matrix(degree, time, query_time, start, M, u_t, dt, time_index_offset))
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

    vector<bspline_trajectory::nbs_pva_state_1d> bspline_trajectory::get_nbspline_1d_all(
        int degree, vector<t_p_sc> time, vector<double> cp, double interval, t_p_sc start)
    {
        vector<bspline_trajectory::nbs_pva_state_1d> v_s;
        // According to the paper only able to calculate to degree 3 which is order 4
        if (degree > 3)
            return v_s;
        if (cp.empty() || time.empty())
            return v_s;
        if (time.size() != cp.size() + (degree-1))
        {
            std::cout << KRED << "time vector size not correct!" << KNRM << std::endl;
            return v_s;
        }

        double total_duration = duration<double>(time.back() - time.front()).count();
        int divisions = (int)floor(total_duration / interval);

        for (int i = 1; i < divisions; i++)
        {
            int k = degree + 1;
            Eigen::MatrixXd M;
            double u_t, dt;
            int time_index_offset;

            t_p_sc query_time = time.front() + milliseconds(i * (int)round(interval*1000));
            
            nbs_pva_state_1d s;

            if (!assemble_M_ut_dt_matrix(degree, time, query_time, start, M, u_t, dt, time_index_offset))
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

    bspline_trajectory::nbs_pva_state_1d bspline_trajectory::get_nbspline_1d_w_prior(
        double dt, int time_index_offset, vector<double> cp, double u_t, Eigen::MatrixXd M, int k)
    {
        bspline_trajectory::nbs_pva_state_1d s;

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

    bspline_trajectory::nbs_pva_state_3d bspline_trajectory::get_nbspline_3d(
        int degree, vector<t_p_sc> time, vector<Eigen::Vector3d> cp, t_p_sc query_time, t_p_sc start)
    {
        nbs_pva_state_3d ss;
        ss.rts = query_time;

        // According to the paper only able to calculate to degree 3 which is order 4
        if (degree > 3)
            return ss;
        if (cp.empty() || time.empty())
        {
            std::cout << KRED << "cp and time vector empty!" << KNRM << std::endl;
            return ss;
        }
        if (time.size() != cp.size() + (degree-1))
        {
            std::cout << KRED << "time vector size not correct!" << KNRM << std::endl;
            return ss;
        }
        int k = degree + 1;
        Eigen::MatrixXd M;
        double u_t, dt;
        int time_index_offset;

        if (!assemble_M_ut_dt_matrix(degree, time, query_time, start, M, u_t, dt, time_index_offset))
        {
            std::cout << KRED << "Cannot assemble time du matrix!" << KNRM << std::endl;
            return ss;
        }

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


        bspline_trajectory::nbs_pva_state_1d x = get_nbspline_1d_w_prior(
            dt, time_index_offset, rv.xcp, u_t, M, k);
        bspline_trajectory::nbs_pva_state_1d y = get_nbspline_1d_w_prior(
            dt, time_index_offset, rv.ycp, u_t, M, k);
        bspline_trajectory::nbs_pva_state_1d z = get_nbspline_1d_w_prior(
            dt, time_index_offset, rv.zcp, u_t, M, k);

        ss.pos = Eigen::Vector3d(x.pos, y.pos, z.pos);
        ss.vel = Eigen::Vector3d(x.vel, y.vel, z.vel);
        ss.acc = Eigen::Vector3d(x.acc, y.acc, z.acc);

        return ss;
    }

    vector<bspline_trajectory::nbs_pva_state_3d> bspline_trajectory::get_nbspline_3d_all(
        int degree, vector<t_p_sc> time, vector<Eigen::Vector3d> cp, double interval, t_p_sc start)
    {
        vector<bspline_trajectory::nbs_pva_state_3d> v_ss;

        // According to the paper only able to calculate to degree 3 which is order 4
        if (degree > 3)
            return v_ss;
        if (cp.empty() || time.empty())
            return v_ss;
        if (time.size() != cp.size() + (degree-1))
        {
            std::cout << KRED << "time vector size not correct!" << KNRM << std::endl;
            return v_ss;
        }
        int k = degree + 1;

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
            bspline_trajectory::nbs_pva_state_3d ss;
            Eigen::MatrixXd M;
            double u_t, dt;
            int time_index_offset;

            t_p_sc query_time = time.front() + milliseconds(i * (int)round(interval*1000));

            if (!assemble_M_ut_dt_matrix(degree, time, query_time, start, M, u_t, dt, time_index_offset))
                return v_ss;
            
            ss.rts = query_time;
            bspline_trajectory::nbs_pva_state_1d x = get_nbspline_1d_w_prior(
                dt, time_index_offset, rv.xcp, u_t, M, k);
            bspline_trajectory::nbs_pva_state_1d y = get_nbspline_1d_w_prior(
                dt, time_index_offset, rv.ycp, u_t, M, k);
            bspline_trajectory::nbs_pva_state_1d z = get_nbspline_1d_w_prior(
                dt, time_index_offset, rv.zcp, u_t, M, k);

            ss.pos = Eigen::Vector3d(x.pos, y.pos, z.pos);
            ss.vel = Eigen::Vector3d(x.vel, y.vel, z.vel);
            ss.acc = Eigen::Vector3d(x.acc, y.acc, z.acc);

            v_ss.push_back(ss);
        }

        return v_ss;
    }

    void bspline_trajectory::assemble_u_p_matrix(
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

    bool bspline_trajectory::assemble_M_ut_dt_matrix(
        int d, vector<t_p_sc> t, t_p_sc q, t_p_sc s, 
        Eigen::MatrixXd& M, double& u_t, double& dt, int& t_i_o)
    {
        int k = d + 1;
        std::pair<t_p_sc, t_p_sc> t_i;
        if (!check_query_time(d, t, q, t_i, t_i_o))
        {
            std::cout << KRED << "query not inside time vector = " << duration<double>(q - t[t.size()-1]).count() << KNRM << std::endl;
            return false;
        }

        vector<double> t_t;
        // std::cout << "time_trim vector";
        for (int i = 0; i < k + (d - 1); i++)
        {
            double specific_time = duration<double>(t[t_i_o + i] - s).count();
            t_t.push_back(specific_time);
            // std::cout << " " << specific_time;
        }
        // std::cout << std::endl;

        M = create_general_m(d, t_t);

        // u_t = (query_time - t_i.first) / (t_i.second - t_i.first)
        double numerator = duration<double>(q - t_i.first).count();
        double denominator = duration<double>(t_i.second - t_i.first).count();
        // Only considering [0,1) hence not including 1
        u_t = numerator / denominator;
        dt = denominator;

        return true;
    }
}