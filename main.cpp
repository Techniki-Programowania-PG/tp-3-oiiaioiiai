#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <matplot/matplot.h>
#include <complex>
#include <pybind11/complex.h>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace matplot;

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("sinus", [](double A, double f, double t0, double t1, int l) {


        std::vector<double> x = linspace(t0, t1, l);
        std::vector<double> y;

        for (int i = 0;i < x.size();i++) {
            y.push_back(std::sin(f*x[i]) * A);
        }

        plot(x, y);
        xlabel("os x");
        ylabel("os y");

        grid(on);
        show();
        save("sinus", "png");


        return y;
    });

    m.def("cosinus", [](double A, double f, double t0, double t1, int l) {


        std::vector<double> x = linspace(t0, t1, l);
        std::vector<double> y;

        for (int i = 0;i < x.size();i++) {
            y.push_back(std::cos(f*x[i]) * A);
        }

        plot(x, y);
        xlabel("os x");
        ylabel("os y");

        grid(on);
        show();

        save("cosinus", "png");
    });

    m.def("prostokatna", [](double A, double f, double t0, double t1, int l) {

        std::vector<double> x = linspace(t0, t1, l);
        std::vector<double> y;

        for (int i = 0;i < x.size();i++) {
            if((std::cos(f * x[i]) * A) >= 0){
                y.push_back(A);
            }
            else {
                y.push_back(-A);
            }
        }

        plot(x, y);
        xlabel("os x");
        ylabel("os y");

        grid(on);
        show();

        save("prostokatna", "png");
    });


    m.def("pila", [](double A, double f, double t0, double t1, int l) {

        std::vector<double> x = linspace(t0, t1, l);
        std::vector<double> y;

        for (int i = 0;i < x.size();i++) {
            y.push_back((fmod((A * x[i] * f) , 2))-1);
        }

        plot(x, y);
        xlabel("os x");
        ylabel("os y");

        grid(on);
        show();

        save("piloksztaltna", "png");
    });

    m.def("Fourier", [](std::vector<double> wejscie, double t0, double t1, int l) {

        std::vector<double> x = linspace(t0, t1, l);
        std::vector<std::complex<double>> y;
        std::vector<double> modul;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> zespolona;
        double niei;

        for (int a = 0;a < l;a++) {
            zespolona = 0;
            for (int b = 0;b < l;b++) {
                niei = 2.0 * pi * double(a) * double(b);
                zespolona += wejscie[b] * std::exp((- i * niei) / double(l));
            }
            y.push_back(zespolona);
            modul.push_back(std::abs(zespolona));
        }
        

        
        plot(x, modul);
        xlabel("os x");
        ylabel("os y");

        grid(on);
        show();

        save("fourier", "png");

        return y;
    });

    m.def("inverse", [](std::vector<std::complex<double>> wejscie, double t0, double t1, int l) {

        std::vector<double> x = linspace(t0, t1, l);
        std::vector<std::complex<double>> y;
        std::vector<double> modul;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> zespolona;
        double niei;

        for (int a = 0;a < l;a++) {
            zespolona = 0;
            for (int b = 0;b < l;b++) {
                niei = 2.0 * pi * double(a) * double(b);
                zespolona += wejscie[b] * std::exp((i * niei) / double(l));
            }
            zespolona = zespolona / double(l);
            y.push_back(zespolona);
            modul.push_back(std::real(zespolona));
        }



        plot(x, modul);
        xlabel("os x");
        ylabel("os y");

        grid(on);
        show();

        save("inverse", "png");
    });

    m.def("jeden", [](std::vector<double> wejscie,double t0, double t1, int l) {

        std::vector<double> x = linspace(t0, t1, l);
        std::vector<double> y;
        std::vector<double> filtr;

        int iks=1;
        for (int i = 0;i < wejscie.size();i++) {
            filtr.push_back(iks*(pow(-1,i+1)));
            iks += 1;
        }
        
        double suma;
        for (int i = 0;i < wejscie.size();i++) {
            suma = 0.0;
            for (int j = 0;j < i + 1;j++) {
                suma += wejscie[j]*filtr[i-j];
            }
            y.push_back(suma);
        }


        plot(x,y);
        xlabel("os x");
        ylabel("os y");

        grid(on);
        show();

        save("filtracja 1d", "png");
    });

    m.def("dwa", [](double t0, double t1) {
        t1 += 1;
        int wielkosc = t1+1;

        std::vector<double> x = linspace(t0, t1, wielkosc);
        std::vector<double> y = linspace(t0, t1, wielkosc);

        std::vector<std::vector<double>> X(wielkosc, std::vector<double>(wielkosc));
        std::vector<std::vector<double>> Y(wielkosc, std::vector<double>(wielkosc));

        for (int i = 0; i < wielkosc -1; ++i) {
            for (int j = 0; j < wielkosc -1; ++j) {
                X[i][j] = x[j];
                Y[i][j] = y[i];
            }
        }

        std::vector<std::vector<double>> dane(wielkosc, std::vector<double>(wielkosc));
        double srodek = (wielkosc - 1) / 2.0;
        double szerokosc = wielkosc / 5.0;

        for (int i = 0; i < wielkosc; i++) {
            for (int j = 0; j < wielkosc; j++) {
                double dx = i - srodek;
                double dy = j - srodek;
                dane[i][j] = std::exp(-(dx * dx + dy * dy) / (2 * szerokosc * szerokosc));
            }
        }

        std::vector<std::vector<double>> filtr{
            {0, -1, 0},
            {-1, 5, -1},
            {0, -1, 0}
        };

        std::vector<std::vector<double>> Z(wielkosc, std::vector<double>(wielkosc, 0.0));

        for (int i = 1; i < wielkosc -1; i++) {
            for (int j = 1; j < wielkosc -1; j++) {
                double suma = 0.0;
                for (int k = -1; k <= 1; k++) {
                    for (int m = -1; m <= 1; m++) {
                        suma += filtr[k + 1][m + 1] * dane[i + k][j + m];
                    }
                }
                Z[i][j] = suma;
            }
        }

        surf(X, Y, Z);
        xlabel("os x");
        ylabel("os y");
        zlabel("os z");
        colorbar();
        grid(on);
        show();

        save("filtracja 2d", "png");

        return Z;
    });


    m.attr("__version__") = "dev";
}
