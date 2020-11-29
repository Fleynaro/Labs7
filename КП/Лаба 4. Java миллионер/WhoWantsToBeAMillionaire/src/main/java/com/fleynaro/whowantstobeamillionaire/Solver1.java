/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.fleynaro.whowantstobeamillionaire;

/**
 *
 * @author user
 */
import static java.lang.Math.*;

interface Functions{
    double value(double x);
}
interface FunctionsXYZ{
    double value(double x, double y, double z);
}

public class Solver1 {


    // Ф-ции наших вариантов
    public static Functions function_A = (x) ->  50*(x + 1); //30*(-x + 0.5);
    public static Functions function_B = (x) ->  -x*x +2; //x*x + 1;
    public static Functions function_C = (x) ->  2*x+1; //x + 2;

    // Упр, У'пр, У"пр
    public static Functions function_Y = (x) -> 1 + x + 10*log(18 + 1)*pow(x, 3)*pow((1 - x), 3);
    public static Functions dFunction_Y = (x) -> 1 - 30*log(18 + 1)*pow(x, 2)*pow((x - 1), 2)*(2*x - 1);
    public static Functions ddFunction_Y = (x) -> -60*log(18 + 1)*x*(x - 1)*(5*pow(x, 2) - 5*x + 1);

    //F(x) = ..., Y"(x) = f = F(x) - A(x)*z - B(x)*y + C(x)*sin(y)
    public static Functions function_F = (x) -> ddFunction_Y.value(x) + function_A.value(x)*dFunction_Y.value(x) -
            function_B.value(x)*function_Y.value(x) + function_C.value(x)*sin(function_Y.value(x));

    public static FunctionsXYZ function_f = (x, y, z) -> function_F.value(x) - function_A.value(x)*z +
            function_B.value(x)*y - function_C.value(x)*sin(y);

    public static double Z0;
    public static double delta;
    public static double alpha0 = 0;
    public static double alpha1 = 3;

    public static double[] stepValue(double x, double y, double z, double h) {

        double K11;
        double K12;
        double K13;
        double K14;
        double K21;
        double K22;
        double K23;
        double K24;

        K11 = h * z;
        K21 = h * function_f.value(x, y, z);
        K12 = h * (z + K21/2.0);
        K22 = h * function_f.value(x + h/2.0, y + K11/2.0, z + K21/2.0);
        K13 = h * (z + K22/2.0);
        K23 = h * function_f.value(x + h/2.0, y + K12/2.0, z + K22/2.0);
        K14 = h * (z + K23);
        K24 = h * function_f.value(x + h, y + K13, z + K23);

        y = y + (K11 + 2*K12 + 2*K13 + K14)/6.0;
        z = z + (K21 + 2*K22 + 2*K23 + K24)/6.0;

        return new double[] {y, z};
    }

    public static double RKMethodPrint() {

        delta = 0;
        double x = 0;
        double y = 1;
        double _y = 0;
        double z = Z0;
        double _z;

        System.out.println();
        System.out.println("|x      |y(x)   |Yпр(x) |z(x)     |delta       |");

        while (x < 1)
        {
            double h = 0.1;
            double y1;
            double y2;
            double y3;
            double z1;

            do {
                _y = y;
                _z = z;
                z1 = z;
                h = h/2.0;


                double[] yz1 = stepValue(x, _y, z1, h);
                y1 = yz1[0];
                double[] yz2 = stepValue(x, _y, _z, h/2.0);
                _z = yz2[1]; y2 = yz2[0];
                double[] yz3 = stepValue(x, y2, _z, h/2.0);
                y3 = yz3[0];

            } while (abs(y1 - y3) > 0.00001);

            System.out.printf("|%7.5f|%7.5f|%7.5f|%9.6f|%9e|", x, y, function_Y.value(x), z, abs(function_Y.value(x) - y));
            System.out.println();

            _y = y;

            double[] yz = stepValue(x, y, z, h);
            y = yz[0]; z = yz[1];

            x = x + h;
        }
        return _y;
    }


    public static double RKMethod() {

        delta = 0;
        double err;
        double x = 0;
        double y = 1;
        double _y = 0;
        double z = Z0;
        double _z;

        while (x < 1)
        {
            double h = 0.1;
            double y1;
            double y2;
            double y3;
            double z1;

            do {

                h = h/2.0;
                _y = y;
                _z = z;
                z1 = z;

                double[] yz1 = stepValue(x, _y, z1, h);
                y1 = yz1[0];
                double[] yz2 = stepValue(x, _y, _z, h/2.0);
                y2 = yz2[0]; _z = yz2[1];
                double[] yz3 = stepValue(x, y2, _z, h/2.0);
                y3 = yz3[0];

            } while (abs(y1 - y3) > 0.00001);

            _y = y;

            err = abs(function_Y.value(x) - y);
            if (err > delta)
                delta = err;

            double[] yz = stepValue(x, y, z, h);
            y = yz[0]; z = yz[1];

            x = x + h;
        }
        return _y;
    }


    public static void shootMethod() {


        int iteration = 0;
        double y;

        System.out.println("Метод стрельб");
        System.out.println("|Itr |z(0)   |y(1)    |Delta       |");

        do {

            Z0 = (alpha1 + alpha0) / 2.0;
            y = RKMethod();

            if (y > 2) {
                alpha1 = Z0;
            } else{
                alpha0 = Z0;
            }

            iteration += 1;

            System.out.printf("|%4d|%7.5f|%8.6f|%12.10f|", iteration, Z0, y, delta);
            System.out.println();


        } while (abs(y - 2) > 0.0001);
        RKMethodPrint();
    }
}
