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

namespace plt = matplotlibcpp;

int main()
{
    std::random_device dev;
    std::mt19937 generator(dev());

    double maximum_variation = 2.0;
    std::uniform_real_distribution<double> dis(1.0, 3.0);
    std::uniform_real_distribution<double> dis_3d(-maximum_variation, maximum_variation);
    std::uniform_real_distribution<double> dis_3d_h(1.0, 2.0);

    bspline_trajectory nb;

    int order = 3;

    /** @brief Key components cp_size and time_point_size relationship **/
    int cp_size = 20;
    int time_point_size = cp_size + order - 1;

    if (time_point_size < 0)
        return -1;

    /** @brief Creation of the knot vector **/
    vector<double> t; // time vector 
    vector<Eigen::Vector3d> cp_3d; // control points 3d

    double t_s = 0.0; // total accumulated time
    std::cout << "time_vector =";
    t.push_back(t_s);
    std::cout << " " << t_s;
    for (int i = 1; i < time_point_size; i++)
    {
        t_s += dis(generator);
        // t_s += 1.0;
        t.push_back(t_s);
        std::cout << " " << t_s;
    }
    std::cout << std::endl;


    /** @brief Creation of the control point vector **/
    Eigen::Vector3d s_cp = 
        Eigen::Vector3d(dis_3d(generator), dis_3d(generator), dis_3d_h(generator));
    Eigen::Vector3d e_cp = 
        Eigen::Vector3d(dis_3d(generator), dis_3d(generator), dis_3d_h(generator));
    std::cout << "control point vector" << std::endl;
    for (int i = 0; i < order; i++)
    {
        cp_3d.push_back(s_cp);
        std::cout << s_cp.transpose() << std::endl;
    }
    for (int i = 0; i < cp_size - 2*order; i++)
    {
        Eigen::Vector3d rand_value = 
            Eigen::Vector3d(dis_3d(generator), dis_3d(generator), dis_3d_h(generator));
        std::cout << rand_value.transpose() << std::endl;
        cp_3d.push_back(rand_value);
    }
    for (int i = 0; i < order; i++)
    {
        cp_3d.push_back(e_cp);
        std::cout << e_cp.transpose() << std::endl;
    }

    std::uniform_real_distribution<double> dis_t(t.front(), ceil(t.back()));
    double q_t = dis_t(generator); // query time

    std::pair<double,double> t_i;

    /** @brief testing of check_query_time, find the appropriate time point pair for evaluation**/
    t_s = t[(cp_size-1)-(order-1)];
    std::cout << "total nbspline time: " << KYEL << t_s << "s" << KNRM << std::endl;
    std::pair<vector<double>,vector<Eigen::Vector3d>> three_d_pos_time; // 1d position vector with time
    time_point<std::chrono::system_clock> t_s_nb = system_clock::now();
    while (duration<double>(system_clock::now() - t_s_nb).count() < t_s)
    {
        time_point<std::chrono::system_clock> t_s_gn3 = system_clock::now();
        q_t = duration<double>(system_clock::now() - t_s_nb).count();
        bspline_trajectory::nbs_pva_state_3d state_3d;
        state_3d = nb.get_nbspline_3d(order, t, cp_3d, q_t);
        three_d_pos_time.second.push_back(state_3d.pos);
        three_d_pos_time.first.push_back(q_t);
        auto t_gn3 = duration<double>(system_clock::now() - t_s_gn3).count() * 1000;
        std::cout << "[" << KYEL << q_t << KNRM << 
            "] get_nbspline_3d " << KGRN << t_gn3 << "ms" << KNRM << std::endl;

        sleep_for(milliseconds(100));
    }

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
    plt::named_plot("x/m", three_d_pos_time.first, rv.xcp, "b--");
    plt::named_plot("y/m", three_d_pos_time.first, rv.ycp, "r--");
    plt::named_plot("z/m", three_d_pos_time.first, rv.zcp, "y--");

    // Just for visualization
    vector<double> t_trim;
    for(int i = 0; i <= (cp_size-1)-(order-1); i++)
        t_trim.push_back(t[i]);

    int plot_segments = 10;
    double steps = (2 * maximum_variation) / (double)(plot_segments-1);
    for(int i = 0; i < (int)t_trim.size(); i++)
    {
        double accumulated = - maximum_variation;
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
    
    string title = std::to_string(order) + " order 3d non-uniform-bspline";
    plt::title(title); // add graph title
    plt::show();

    return 0;
}