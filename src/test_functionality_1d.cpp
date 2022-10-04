/*
* test_functionality_1d.hpp
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
#include "matplotlibcpp.h"
#include <iostream>
#include <random>
#include <random>
#include <thread>

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

using namespace nbspline;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds

typedef time_point<std::chrono::system_clock> t_p_sc; // giving a typename

namespace plt = matplotlibcpp;

int main(int argc, char** argv)
{
    std::random_device dev;
    std::mt19937 generator(dev());

    int fast_forward = 0;
    if (argc == 1)
        std::cout << KRED << "[no input] default no fast forward" << KNRM << std::endl;
    else if (strcmp(argv[1], "fast") == 0)
    {
        fast_forward = 1;
        std::cout << KGRN << "fast-forward activated" << KNRM << std::endl;
    }
    else if (strcmp(argv[1], "faster") == 0)
    {
        fast_forward = 2;
        std::cout << KGRN << "even fast-forward activated" << KNRM << std::endl;
    }
    else
    {
        std::cout << KRED << "[no proper value] default no fast forward" << KNRM << std::endl;
    }

    double maximum_variation = 2.0;
    std::uniform_real_distribution<double> dis(1.0, 3.0);
    std::uniform_real_distribution<double> dis_1d(-maximum_variation, maximum_variation);

    bspline_trajectory nb;

    int run_interval_ms = 100;

    int order = 3;
    int k = order+1;

    /** @brief Key components cp_size and time_point_size relationship **/
    int cp_size = 10;
    int time_point_size = cp_size + (order-1);

    if (time_point_size < 0)
        return -1;

    /** @brief Creation of the knot vector **/
    vector<t_p_sc> t; // time vector in chronos time point
    vector<double> cp_1d; // control points 1d

    /** @brief Get the current start time of the bspline **/
    t_p_sc t_start = system_clock::now() + seconds(1);
    
    double t_s = 0.0; // total accumulated time
    int t_a_ms = 0; // total accumulated time in ms
    std::cout << "time_vector =";
    
    for (int i = 0; i < time_point_size; i++)
    {
        if (i < order)
        {
            t.push_back(t_start);
            std::cout << " " << t_s;
            continue;
        }

        if (i > time_point_size - order)
        {
            t.push_back(t.back());
            std::cout << " " << t_s;
            continue;
        }

        /** @brief Uniform distribution **/
        // double t_single_ms = round(dis(generator) * 1000);
        // int t_ms = (int)t_single_ms; // milliseconds
        // t_s += t_single;
        // t_a_ms += t_ms;
        // t.push_back(t_start + milliseconds(t_a_ms);
        // std::cout << " " << t_s;
        
        /** @brief Non-uniform distribution **/
        double t_single_ms = round(dis(generator) * 1000);
        int t_ms = (int)t_single_ms; // milliseconds
        t_s += t_single_ms / 1000;
        t_a_ms += t_ms;
        t.push_back(t_start + milliseconds(t_a_ms));
        std::cout << " " << t_s;   
    }
    std::cout << std::endl;

    /** @brief Creation of the control point vector **/
    double s_cp = dis_1d(generator);
    double e_cp = dis_1d(generator);
    std::cout << "control point vector =";
    for (int i = 0; i < cp_size; i++)
    {
        if (i < order)
        {
            cp_1d.push_back(s_cp);
            std::cout << " " << s_cp;
            continue;
        }

        if (i >= cp_size - order)
        {
            cp_1d.push_back(e_cp);
            std::cout << " " << e_cp;
            continue;
        }

        double rand_value = dis_1d(generator);
        std::cout << " " << rand_value;
        cp_1d.push_back(rand_value);
    }
    std::cout << std::endl;

    std::pair<t_p_sc, t_p_sc> t_i;
    int time_index_offset = 0;

    /** @brief testing of check_query_time, find the appropriate time point pair for evaluation**/
    t_p_sc t_s_cqt = system_clock::now();
    if (!nb.check_query_time(order, t, t_start + milliseconds(1), t_i, time_index_offset))
    {
        std::cout << "check_query_time fail, query time not inside time vector" << std::endl;
        return -1;
    }
    auto t_cqt = duration<double>(system_clock::now() - t_s_cqt).count() * 1000;
    std::cout << "check_query_time " << KGRN << t_cqt << "ms" << KNRM << std::endl;

    /** @brief testing create_general_m, creation of the M matrix representing the basis of the spline**/
    t_p_sc t_s_cgm = system_clock::now();
    vector<double> time_trim;
    std::cout << "time_trim vector test";
    for (int i = 0; i < k + (order-1); i++)
    {
        double specific_time = duration<double>(t[time_index_offset+i] - t_start).count();
        time_trim.push_back(specific_time);
        std::cout << " " << specific_time;
    }
    std::cout << std::endl;
    Eigen::MatrixXd m = nb.create_general_m(order, time_trim);
    auto t_cgm = duration<double>(system_clock::now() - t_s_cgm).count();
    std::cout << "create_general_m " << KGRN << t_cgm * 1000 << "ms" << KNRM << std::endl;
    std::cout << m << std::endl;


    auto us_interval = std::chrono::duration_cast<std::chrono::milliseconds>(t_start - system_clock::now());
    sleep_for(us_interval + milliseconds(1));

    std::cout << "total nbspline time: " << KYEL << t_s << "s" << KNRM << std::endl;
    std::pair<vector<double>,vector<double>> one_d_pos_time; // 1d position vector with time
    vector<double> one_d_vel, one_d_acc;

    t_p_sc now = system_clock::now();
    t_p_sc t_s_gnt = now;
    if (fast_forward != 2)
    {
        while (duration<double>(now - t_start).count() < t_s)
        {
            t_p_sc t_s_gn1 = system_clock::now();
            
            bspline_trajectory::nbs_pva_state_1d state_1d;
            state_1d = nb.get_nbspline_1d(order, t, cp_1d, now, t_start);
            one_d_pos_time.second.push_back(state_1d.pos);

            one_d_pos_time.first.push_back(duration<double>(now - t_start).count());
            one_d_vel.push_back(state_1d.vel);
            one_d_acc.push_back(state_1d.acc);

            auto t_gn1 = duration<double>(system_clock::now() - t_s_gn1).count();
            std::cout << "[" << KYEL << duration<double>(now - t_start).count() << KNRM << 
                "] get_nbspline_1d " << KGRN << t_gn1 * 1000 << "ms" << KNRM << std::endl;

            if (fast_forward != 1)
            {
                sleep_for(milliseconds(run_interval_ms));
                now = system_clock::now();
            }
            else
                now += milliseconds(run_interval_ms);

        }
    }
    else
    {
        vector<bspline_trajectory::nbs_pva_state_1d> state_1d_vector;
        state_1d_vector = nb.get_nbspline_1d_all(order, t, cp_1d, run_interval_ms/1000.0, t_start);
        for (auto state_1d : state_1d_vector)
        {
            one_d_pos_time.second.push_back(state_1d.pos);
            one_d_pos_time.first.push_back(duration<double>(state_1d.rts - t_start).count());
            one_d_vel.push_back(state_1d.vel);
            one_d_acc.push_back(state_1d.acc);
        }
    }

    std::cout << "time taken for 1d bspline size (" << one_d_pos_time.first.size() << ") calculation is " << 
        duration<double>(system_clock::now() - t_s_gnt).count() * 1000 << "ms" << std::endl;

    // vector<double> velocity_vect, acceleration_vect;
    // velocity_vect.push_back(0.0);
    // for (int i = 1; i < (int)one_d_pos_time.second.size(); i++)
    // {
    //     velocity_vect.push_back(
    //         (one_d_pos_time.second[i]-one_d_pos_time.second[i-1])/((double)run_interval_ms/1000.0));
    //     if (i >= 2)
    //         acceleration_vect.push_back(
    //             (velocity_vect[i-1]-velocity_vect[i-2])/((double)run_interval_ms/1000.0));
    // }
    // acceleration_vect.push_back(acceleration_vect[acceleration_vect.size()-1]);
    // acceleration_vect.push_back(acceleration_vect[acceleration_vect.size()-1]);

    // Set the size of output image to 1200x780 pixels
    plt::figure_size(980, 460);
    // plot a red dashed line from given x and y data.
    plt::named_plot("pos", one_d_pos_time.first, one_d_pos_time.second, "b*");
    plt::named_plot("vel", one_d_pos_time.first, one_d_vel, "r*");
    plt::named_plot("acc", one_d_pos_time.first, one_d_acc, "y*");
    // plt::named_plot("check_velocity", one_d_pos_time.first, velocity_vect, "c--");
    // plt::named_plot("check_acceleration", one_d_pos_time.first, acceleration_vect, "g--");

    // Just for visualization
    vector<double> t_trim;
    for(int i = order-1; i < time_point_size - (order-1); i++)
        t_trim.push_back(duration<double>(t[i] - t_start).count());

    int plot_segments = 10;
    double start_range = 1.5 * maximum_variation;
    double steps = (2 * start_range) / (double)(plot_segments-1);
    for(int i = 0; i < (int)t_trim.size(); i++)
    {
        double accumulated = - start_range;
        vector<double> t_plot(plot_segments, t_trim[i]);
        vector<double> l_plot;
        for (int j = 0; j < plot_segments; j++)
        {
            l_plot.push_back(accumulated);
            accumulated += steps;
        }
        plt::plot(t_plot, l_plot, "m--");
    }

    // Enable legend.
    plt::legend();
    
    string title = std::to_string(order) + " order 1d non-uniform-bspline";
    plt::title(title); // add graph title
    plt::show();

    return 0;
}