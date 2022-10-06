# Documentation for B-splines

## Some resources
- https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html
- https://xiaoxingchen.github.io/2020/03/02/bspline_in_so3/general_matrix_representation_for_bsplines.pdf
- https://arxiv.org/pdf/1703.01416.pdf

---

## Introduction
B-spline or Basis Spline, is a piecewise polynomial function of degree in a variable. The variables for B-splines are in the form of `knots` and `control points`.

We can define the knot vector $U$ to be an ascending vector of `time points` where in this implementation, we use `time_point<std::chrono::system_clock>`. (etc 0, 0.5, 1.0, 1.5, 2.0)

`Control points` in the simplest form would be a 1d points scattered across, hence the plot for B-spline would be in `x(m) - t(s)` domain

$$
P(t) = \sum_{i=0}^n N_{i,k}(t) \ p_i \ , \\ 
s.t \ \ \ t_{k-1} \le t \le t_{n+1}
$$

Basis function is represented by recurrence relationship (see de boor algorithm below)

Basis function: 
$$
N_{i,k} = \frac{t - t_i}{t_{i+k-1} - t_i} N_{i,k-1}(t) + \frac{t_{i+k} - t}{t_{i+k} - t_{i+1}} N_{i+1,k-1}(t)
$$
$k$ represent the `order` (`degree` = `order-1`), $i$ represents the `control point` being used in evaluating the spline

For example, i = 0, k = 2, then the basis function is:
$$
N_{0,2} = \frac{t - t_0}{t_{1} - t_0} N_{0,1}(t) + \frac{t_{2} - t}{t_{2} - t_{1}} N_{1,1}(t)
$$

In piecewise form, 

| piecewise form | bspline variables |
| :-: | :-: |
| [<img src="media/bspline_pp.png" width="1000"/>](media/bspline_pp.png) | [<img src="media/bsplinebasics.png" width="600"/>](media/bspline_pp.png) |

A good examples of bsplines would be https://www.ibiblio.org/e-notes/Splines/basis.html

---
## Cox-de Boor Algorithm
It is a fast and numerically stable way for finding a point on a B-spline curve given a u in the domain, follow the notes from https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/de-Boor.html if there are things you want to clarify



We can take a look at `MTU`'s course on Bsplines, the example given is a 3rd degree B-spline curve defined by 
- 7 control points `p0, p1 ... p6` 
- 11 knots `u0 to u3 = 0, u4 = 0.25, u5 = 0.5, u6 = 0.75, u7 to u10 = 1`

And the example here is that we want to find `p(0.4)`, note that it lies between knots `u4` and `u5`

The notation here for the `de-Boor algorithm` is that in $P_{i,j}$ (points from the basis curve)
- $j$ represents the columns (we have to reduce the columns to 1 to find the final point)
- $i$ represents the points (after being altered when the degree is raised)

| Correct Formulation | Wrong example (not using formulation) |
| :-: | :-: |
| [<img src="media/de-boor-algo.jpg" width="300"/>](media/de-boor-algo.jpg) | [<img src="media/de-boor-ex-1.jpg" width="300"/>](media/de-boor-ex-1.jpg) |

**Wrong example** shown in https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/de-Boor.html this does not follow de-boor algorithm


We would start with control points $p_{1,0}, p_{2,0}, p_{3,0}, p_{4,0}$ and this would be the same as their basis function counterpart.
$$
N_{1,0} = p_{1,0}
$$

### From column 1 to 2:

~~We can find the points for `column 2` in the **figure above** using the equation below~~

Do not use this example below

$$
p_{4,1} = (1 - u_{4,1}) \ p_{3,0} + u_{4,1} \ p_{4,0} \\
p_{3,1} = (1 - u_{3,1}) \ p_{2,0} + u_{3,1} \ p_{3,0} \\
p_{2,1} = (1 - u_{2,1}) \ p_{1,0} + u_{2,1} \ p_{2,0}
$$

$$
u_{4,1} = \frac{t - t_4}{t_{4+3} - t_4} = 0.2 \\
u_{3,1} = \frac{t - t_3}{t_{3+3} - t_3} = 0.53 \\
u_{2,1} = \frac{t - t_2}{t_{2+3} - t_2} = 0.8 
$$

### From column 2 to 3:

~~We can find the points for `column 3` in the **figure above** using the equation below~~

Do not use this example below

$$
p_{4,2} = (1 - a_{4,2}) \ p_{3,1} + a_{4,2} \ p_{4,1} \\
p_{3,2} = (1 - a_{3,2}) \ p_{2,1} + a_{3,2} \ p_{3,1}
$$

$$
a_{4,2} = \frac{u - u_4}{u_{4+3-1} - u_4} = 0.3 \\
a_{3,2} = \frac{u - u_3}{u_{3+3-1} - u_3} = 0.8
$$

### From column 3 to 4:

~~We can find the points for `column 4` in the **figure above** using the equation below~~

Do not use this example below

$$
p_{4,3} = (1 - u_{4,3}) \ p_{3,2} + u_{4,3} \ p_{4,2}
$$

$$
u_{4,3} = \frac{t - t_4}{t_{4+3-2} - t_4} = 0.6
$$


### Summary
In short de-Boor algorithm substitute the weights from each piecewise curve using a recursive formulation to derive the point in time.

---

## General Representation

**Correction** of derivative equations in `Kaihuai Qin`'s implementation by `Vladyslav Usenko`

- **Position** is expressed as $$c_{k-1} = U^k M^k V^k$$
- `Original` **Velocity** is expressed as $$\frac{d}{du}c_{k-1} = \frac{dU^k}{du} M^k V^k$$
- `Corrected` **Velocity** is expressed as $$\frac{d}{du}c_{k-1} = \frac{1}{(t_{i+1} - t_i)} * \frac{dU^k}{du} M^k V^k$$
- `Original` **Acceleration** is expressed as $$\frac{d^2}{du^2}c_{k-1} = \frac{d^2U^k}{du^2} M^k V^k$$
- `Corrected` **Acceleration** is expressed as $$\frac{d^2}{du^2}c_{k-1} = (\frac{1}{(t_{i+1} - t_i)})^2 \frac{d^2U^k}{du^2} M^k V^k$$

---

| **Wrong** without factor | **Wrong** 3d plot for position |
| :-: | :-: |
|[<img src="media/bspline_wrong.png" width="600"/>](media/bspline_wrong.png.png)|[<img src="media/bspline_3d_wrong.png" width="600"/>](media/bspline_3d_wrong.png)|

| **Correct** pos, vel and acc | **Correct** 3d plot for position |
| :-: | :-: |
|[<img src="media/3rd_degree_1d_nbspline.png" width="600"/>](media/3rd_degree_1d_nbspline.png)|[<img src="media/3rd_degree_3d_nbspline.png" width="600"/>](media/3rd_degree_3d_nbspline.png)|

---