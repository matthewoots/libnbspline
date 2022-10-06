/*
* test_functionality_3d.hpp
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
    std::uniform_real_distribution<double> dis_3d(-maximum_variation, maximum_variation);
    std::uniform_real_distribution<double> dis_3d_h(1.0, 2.0);

    bspline_trajectory nb;

    int run_interval_ms = 100;

    int degree = 3;

    /** @brief Key components cp_size and time_point_size relationship **/
    int cp_size = 10;
    int time_point_size = cp_size + (degree-1);

    if (time_point_size < 0)
        return -1;

    /** @brief Creation of the knot vector **/
    vector<t_p_sc> t; // time vector in chronos time point
    vector<Eigen::Vector3d> cp_3d; // control points 1d

    /** @brief Get the current start time of the bspline **/
    t_p_sc t_start = system_clock::now() + seconds(1);
    
    double t_s = 0.0; // total accumulated time
    int t_a_ms = 0; // total accumulated time in ms
    std::cout << "time_vector =";
    
    for (int i = 0; i < time_point_size; i++)
    {
        if (i < degree)
        {
            t.push_back(t_start);
            std::cout << " " << t_s;
            continue;
        }

        if (i > time_point_size - degree)
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
    Eigen::Vector3d s_cp = 
        Eigen::Vector3d(dis_3d(generator), dis_3d(generator), dis_3d_h(generator));
    Eigen::Vector3d e_cp = 
        Eigen::Vector3d(dis_3d(generator), dis_3d(generator), dis_3d_h(generator));
    std::cout << "control point vector" << std::endl;
    for (int i = 0; i < cp_size; i++)
    {
        if (i < degree)
        {
            cp_3d.push_back(s_cp);
            std::cout << s_cp.transpose() << std::endl;
            continue;
        }

        if (i >= cp_size - degree)
        {
            cp_3d.push_back(e_cp);
            std::cout << e_cp.transpose() << std::endl;
            continue;
        }

        Eigen::Vector3d rand_value = 
            Eigen::Vector3d(dis_3d(generator), dis_3d(generator), dis_3d_h(generator));
        std::cout << rand_value.transpose() << std::endl;
        cp_3d.push_back(rand_value);
    }

    auto us_interval = std::chrono::duration_cast<std::chrono::milliseconds>(t_start - system_clock::now());
    sleep_for(us_interval + milliseconds(1));

    std::cout << "total nbspline time: " << KYEL << t_s << "s" << KNRM << std::endl;
    std::pair<vector<double>,vector<Eigen::Vector3d>> three_d_pos_time; // 1d position vector with time
    
    t_p_sc now = system_clock::now();
    t_p_sc t_s_gnt = now;
    if (fast_forward != 2)
    {
        while (duration<double>(now - t_start).count() < t_s)
        {
            t_p_sc t_s_gn3 = system_clock::now();

            bspline_trajectory::nbs_pva_state_3d state_3d;
            state_3d = nb.get_nbspline_3d(degree, t, cp_3d, now, t_start);
            three_d_pos_time.second.push_back(state_3d.pos);

            three_d_pos_time.first.push_back(duration<double>(now - t_start).count());
            
            auto t_gn3 = duration<double>(system_clock::now() - t_s_gn3).count();
            std::cout << "[" << KYEL << duration<double>(now - t_start).count() << KNRM << 
                "] get_nbspline_3d " << KGRN << t_gn3 * 1000 << "ms" << KNRM << std::endl;

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
        vector<bspline_trajectory::nbs_pva_state_3d> state_3d_vector;
        state_3d_vector = nb.get_nbspline_3d_all(degree, t, cp_3d, run_interval_ms/1000.0, t_start);
        for (auto state_3d : state_3d_vector)
        {
            three_d_pos_time.second.push_back(state_3d.pos);
            three_d_pos_time.first.push_back(duration<double>(state_3d.rts - t_start).count());
        }
    }

    std::cout << "time taken for 3d bspline size (" << three_d_pos_time.first.size() << ") calculation is " << 
        duration<double>(system_clock::now() - t_s_gnt).count() * 1000 << "ms" << std::endl;

    bspline_trajectory::row_vector_3d rv;
    // Reorganize the control points into vectors that can be passed into 1d_bspline function
    for (int i = 0; i < (int)three_d_pos_time.second.size(); i++)
    {
        rv.xcp.push_back(three_d_pos_time.second[i].x());
        rv.ycp.push_back(three_d_pos_time.second[i].y());
        rv.zcp.push_back(three_d_pos_time.second[i].z());
    }

    // Set the size of output image to 1200x780 pixels
    plt::figure_size(980, 460);
    // plot a red dashed line from given x and y data.
    plt::named_plot("x/m", three_d_pos_time.first, rv.xcp, "b*");
    plt::named_plot("y/m", three_d_pos_time.first, rv.ycp, "r*");
    plt::named_plot("z/m", three_d_pos_time.first, rv.zcp, "y*");

    // Just for visualization
    vector<double> t_trim;
    for(int i = degree-1; i < time_point_size - (degree-1); i++)
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
    
    string title = std::to_string(degree) + " degree 3d non-uniform-bspline";
    plt::title(title); // add graph title
    plt::show();

    return 0;
}