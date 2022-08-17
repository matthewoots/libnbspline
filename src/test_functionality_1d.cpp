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

namespace plt = matplotlibcpp;

int main()
{
    std::random_device dev;
    std::mt19937 generator(dev());

    double maximum_variation = 2.0;
    std::uniform_real_distribution<double> dis(1.0, 3.0);
    std::uniform_real_distribution<double> dis_1d(-maximum_variation, maximum_variation);

    bspline_trajectory nb;

    int order = 3;
    int k = order+1;

    /** @brief Key components cp_size and time_point_size relationship **/
    int cp_size = 10;
    int time_point_size = cp_size + order - 1;

    if (time_point_size < 0)
        return -1;

    /** @brief Creation of the knot vector **/
    vector<double> t; // time vector 
    vector<double> cp_1d; // control points 1d

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
    double s_cp = dis_1d(generator);
    double e_cp = dis_1d(generator);
    std::cout << "control point vector =";
    for (int i = 0; i < order; i++)
    {
        cp_1d.push_back(s_cp);
        std::cout << " " << s_cp;
    }
    for (int i = 0; i < cp_size - 2*order; i++)
    {
        double rand_value = dis_1d(generator);
        std::cout << " " << rand_value;
        cp_1d.push_back(rand_value);
    }
    for (int i = 0; i < order; i++)
    {
        cp_1d.push_back(e_cp);
        std::cout << " " << e_cp;
    }
    std::cout << std::endl;

    std::uniform_real_distribution<double> dis_t(t.front(), ceil(t.back()));
    double q_t = dis_t(generator); // query time

    std::pair<double,double> t_i;
    int time_index_offset = 0;

    /** @brief testing of check_query_time, find the appropriate time point pair for evaluation**/
    time_point<std::chrono::system_clock> t_s_cqt = system_clock::now();
    if (!nb.check_query_time(t, q_t, t_i, time_index_offset))
    {
        std::cout << "check_query_time fail, query time not inside time vector" << std::endl;
        return -1;
    }
    auto t_cqt = duration<double>(system_clock::now() - t_s_cqt).count() * 1000;
    std::cout << "check_query_time " << KGRN << t_cqt << "ms" << KNRM << std::endl;


    if (!nb.check_query_time(t, 0.001, t_i, time_index_offset))
    {
        std::cout << "check_query_time fail, query time not inside time vector" << std::endl;
        return -1;
    }
    /** @brief testing create_general_m, creation of the M matrix representing the basis of the spline**/
    time_point<std::chrono::system_clock> t_s_cgm = system_clock::now();
    vector<double> time_trim;
    std::cout << "time_trim vector test";
    for (int i = 0; i < k + (order-1); i++)
    {
        time_trim.push_back(t[time_index_offset+i]);
        std::cout << " " << t[time_index_offset+i];
    }
    std::cout << std::endl;
    Eigen::MatrixXd m = nb.create_general_m(order, time_trim);
    auto t_cgm = duration<double>(system_clock::now() - t_s_cgm).count() * 1000;
    std::cout << "create_general_m " << KGRN << t_cgm << "ms" << KNRM << std::endl;
    std::cout << m << std::endl;

    t_s = t[(cp_size-1)-(order-1)];
    std::cout << "total nbspline time: " << KYEL << t_s << "s" << KNRM << std::endl;
    std::pair<vector<double>,vector<double>> one_d_pos_time; // 1d position vector with time
    vector<double> one_d_vel, one_d_acc;
    time_point<std::chrono::system_clock> t_s_nb = system_clock::now();
    while (duration<double>(system_clock::now() - t_s_nb).count() < t_s)
    {
        time_point<std::chrono::system_clock> t_s_gn1 = system_clock::now();
        q_t = duration<double>(system_clock::now() - t_s_nb).count();
        bspline_trajectory::nbs_pva_state_1d state_1d;
        state_1d = nb.get_nbspline_1d(order, t, cp_1d, q_t);
        one_d_pos_time.second.push_back(state_1d.pos);
        one_d_pos_time.first.push_back(q_t);
        one_d_vel.push_back(state_1d.vel);
        one_d_acc.push_back(state_1d.acc);
        auto t_gn1 = duration<double>(system_clock::now() - t_s_gn1).count() * 1000;
        std::cout << "[" << KYEL << q_t << KNRM << 
            "] get_nbspline_1d " << KGRN << t_gn1 << "ms" << KNRM << std::endl;

        sleep_for(milliseconds(100));
    }

    // Set the size of output image to 1200x780 pixels
    plt::figure_size(980, 460);
    // plot a red dashed line from given x and y data.
    plt::named_plot("pos", one_d_pos_time.first, one_d_pos_time.second, "b--");
    plt::named_plot("vel", one_d_pos_time.first, one_d_vel, "r--");
    plt::named_plot("acc", one_d_pos_time.first, one_d_acc, "y--");

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
    
    string title = std::to_string(order) + " order 1d non-uniform-bspline";
    plt::title(title); // add graph title
    // plt::show();

    return 0;
}