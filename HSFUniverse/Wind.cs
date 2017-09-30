using System;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using Utilities;
using HSFUniverse.DataFiles.HWM;
namespace HSFUniverse
{
    
    public class HWM
    {
        //////
        //////  Horizontal Wind Model 14
        //////
        //////  AUTHORS
        //////    Douglas Drob(0 to ~450+ km, quite-time)
        //////    John Emmert(disturbance winds, DWM Emmert et al., (2008))
        //////    Geospace Science and Technology Branch
        //////    Space Science Division
        //////    Naval Research Laboratory
        //////    4555 Overlook Ave.
        //////    Washington, DC 20375
        //////
        //////  Point of Contact
        //////   douglas.drob @nrl.navy.mil
        //////
        //////   DATE
        //////    July 8, 2014
        //////
        //////
        //////
        //////================================================================================
        ////// Input arguments:
        //////        iyd - year and day as yyddd
        //////        sec - ut(sec)
        //////        alt - altitude(km)
        //////        glat - geodetic latitude(deg)
        //////        glon - geodetic longitude(deg)
        //////        stl - not used
        //////        f107a - not used
        //////        f107 - not used
        //////        ap - two element array with
        //////             ap[1] = not used
        //////             ap[2] = current 3hr ap index
        //////
        ////// Output argument:
        //////        w[1] = meridional wind(m/sec + northward)
        //////        w[2] = zonal wind(m/sec + eastward)
        //////
        //////================================================================================
        static int nmaxhwm = 0;       // maximum degree hwmqt
        static int omaxhwm = 0;       // maximum order hwmqt
        static int nmaxdwm = 0;       // maximum degree hwmqt
        static int mmaxdwm = 0;       // maximum order hwmqt
        static int nmaxqdc = 0;       // maximum degree of coordinate coversion
        static int mmaxqdc = 0;       // maximum order of coordinate coversion
        static int nmaxgeo = 0;       // maximum of nmaxhwm, nmaxqd
        static int mmaxgeo = 0;       // maximum of omaxhwm, nmaxqd

        static Matrix<double> gpbar, gvbar, gwbar; // alfs for geo coordinates
        static Matrix<double> spbar, svbar, swbar;// alfs MLT calculation

        static double glatalf = -1;
        static Matrix<double> xcoeff;         //Coefficients for x coordinate
        static Matrix<double> ycoeff;       //Coefficients for y coordinate
        static Matrix<double> zcoeff;          //Coefficients for z coordinate
        static Matrix<double> sh;              //Array to hold spherical harmonic fuctions

        int nmax0, mmax0;
        // static normalizational coeffiecents
        static Matrix<double> anm, bnm, dnm;
        static Matrix<double> cm, en;
        static Matrix<double> marr, narr;

        public QWM qwm;
        public DWM dwm;


        public HWM()
        {

            int nmax0, mmax0;

            qwm = new QWM();
            dwm = new DWM(out nmaxdwm, out mmaxdwm);

            nmaxgeo = Math.Max(nmaxhwm, nmaxqdc);
            mmaxgeo = Math.Max(omaxhwm, mmaxqdc);

            nmax0 = Math.Max(nmaxgeo, nmaxdwm);
            mmax0 = Math.Max(mmaxgeo, mmaxdwm);

            ALF(nmax0, mmax0);


        }
        public Vector hwm14(int iyd, double sec, double alt, double glat, double glon, double stl, double f107a, double f107, Vector ap)
        {

            Vector w = new Vector(2);

            w = qwm.hwmqt(iyd, sec, alt, glat, glon, stl, f107a, f107);

            if (ap[2] >= 0.0)
            {
                Vector dw = dwm.dwm07(iyd, sec, alt, glat, glon, ap);
                w = w + dw;
            }


            return w;
        }
        public static void alfbasis(int nmax, int mmax, double theta, out Matrix<double> P, out Matrix<double> V, out Matrix<double> W)
        {
            // -------------------------------------------------------------
            // routine to compute vector spherical harmonic basis functions
            // -------------------------------------------------------------

            int n, m;
            double x, y;
            double p00 = 0.70710678118654746;
            P = new Matrix<double>(nmax, mmax);
            V = new Matrix<double>(nmax, mmax);
            W = new Matrix<double>(nmax, mmax);
            P[1, 1] = p00;
            x = Math.Cos(theta);
            y = Math.Sin(theta);
            for (m = 2; m <= mmax+1; m++)
            {
                W[m, m] = cm[m] * P[m - 1, m - 1];
                P[m, m] = y * en[m] * W[m, m];
                for (n = m + 1; n <= nmax+1; n++)
                {
                    W[n, m] = anm[n, m] * x * W[n - 1, m] - bnm[n, m] * W[n - 2, m];
                    P[n, m] = y * en[n] * W[n, m];
                    V[n, m] = narr[n] * x * W[n, m] - dnm[n, m] * W[n - 1, m];
                    W[n - 2, m] = marr[m] * W[n - 2, m];
                }
                W[nmax, m] = marr[m] * W[nmax, m];
                W[nmax+1, m] = marr[m] * W[nmax+1, m];
                V[m, m] = x * W[m, m];
            }
            P[2,1] = anm[2,1] * x * P[1,1];
            V[2,1] = -P[2,2];
            for (n = 3; n <= nmax+1; n++)
            {
                P[n, 1] = anm[n, 1] * x * P[n - 1, 1] - bnm[n, 1] * P[n - 2, 1];
                V[n, 1] = -P[n, 2];
            }

            //return;

        }
        public void ALF(int nmaxin, int mmaxin)
            {
                // -----------------------------------------------------
                // routine to compute static normalization coeffiecents
                // -----------------------------------------------------

                //integer[4], intent(in) :: nmaxin, mmaxin
                int n, m;   // 64 bits to avoid overflow for [m, n] > 60

                int nmax0 = nmaxin;
                int mmax0 = mmaxin;

                anm = new Matrix<double>(nmax0-1, mmax0-1);
                bnm = new Matrix<double>(nmax0-1, mmax0-1);
                cm = new Vector(mmax0);
                dnm = new Matrix<double>(nmax0-1, mmax0-1);
                en = new Vector(nmax0);
                marr = new Vector(mmax0);
                narr = new Vector(nmax0);


                for (n = 1; n <= nmax0; n++)
                {
                    narr[n+1] = n;
                    en[n+1] = Math.Sqrt((double)(n * (n + 1)));
                    anm[n+1, 1] = Math.Sqrt((double)((2 * n - 1) * (2 * n + 1))) / narr[n+1];
                    bnm[n+1, 1] = Math.Sqrt((double)((2 * n + 1) * (n - 1) * (n - 1)) / (2 * n - 3)) / narr[n+1];
                }
                for (m = 1; m <= mmax0; m++)
                {
                    marr[m+1] = m;
                    cm[m+1] = Math.Sqrt((double)(2 * m + 1) /(2 * m * m * (m + 1)));
                    for (n = m + 1; n <= nmax0; n++)
                    {
                        anm[n+1, m+1] = Math.Sqrt((double)((2 * n - 1) * (2 * n + 1) * (n - 1)) / ((n - m) * (n + m) * (n + 1)));
                        bnm[n+1, m+1] = Math.Sqrt((double)((2 * n + 1) * (n + m - 1) * (n - m - 1) * (n - 2) * (n - 1))
                                    / ((n - m) * (n + m) * (2 * n - 3) * n * (n + 1)));
                        dnm[n+1, m+1] = Math.Sqrt((double)((n - m) * (n + m) * (2 * n + 1) * (n - 1)) / ((2 * n - 1) * (n + 1)));
                    }
                }
            //}
        }

        public class QWM
        {

            int nbf;             // Count of basis terms per model level
            int maxn;            // latitude
            int maxs, maxm, maxl;   // seasonal,stationary,migrating
            int maxo;

            int p;                // B-splines order, p = 4 cubic, p=3 quadratic
            int nlev;             // e.g.Number of B-spline nodes
            int nnode;            // nlev + p

            double alttns;           // Transition 1
            double altsym;           // Transition 2
            double altiso;           // Constant Limit
            Vector e1 = new Vector(5);
            Vector e2 = new Vector(5);
            double H = 60.0;

            Matrix<int> nb;           // total number of basis functions @ level
            Matrix<int> order;       // spectral content @ level
            Vector vnode;         // Vertical Altitude Nodes
            public Matrix<double> mparm;      // Model Parameters
            public Matrix<double> tparm;       // Model Parameters

            Vector previous;
            int priornb = 0;

            Matrix<double> fs, fm, fl;
            public Vector bz, bm;

            Vector zwght;
            int lev;

            int cseason = 0;
            int cwave = 0;
            int ctide = 0;

            bool[] content = Enumerable.Repeat(true, 5).ToArray();      // Season/Waves/Tides
            bool[] component = Enumerable.Repeat(true, 2).ToArray();    // Compute zonal/meridional

            //string qwmdefault = "hwm123114.bin";
            //bool qwminit = true;

            Vector wavefactor = new Vector(new List<double>(new double[] { 1.0, 1.0, 1.0, 1.0 }));
            Vector tidefactor = new Vector(new List<double>(new double[] { 1.0, 1.0, 1.0, 1.0 }));

            //HWM parent;

            public QWM()
            {
                int i, j;

                nbf = HWMData.hwm.nbf;
                maxs = HWMData.hwm.maxs;
                maxm = HWMData.hwm.maxm;
                maxl = HWMData.hwm.maxl;
                maxn = HWMData.hwm.maxn;
                nlev = HWMData.hwm.nlev;
                p = HWMData.hwm.p;
                nnode = nlev + p;
                vnode = HWMData.hwm.vnode;
                mparm = HWMData.hwm.mparm;
                e1 = HWMData.hwm.e1;
                e2 = HWMData.hwm.e2;
                nb = HWMData.hwm.nb;
                order = HWMData.hwm.order;
                
                // Calculate the parity relationship permutations

                tparm = new Matrix<double>(nbf, nlev + 1);
                for (i = 1; i < nlev - p +1; i++)
                {
                    int[] orderKarray = new int[order.NumRows];
                    for (int k = 1; k <= order.NumRows; k++)
                    {
                        orderKarray[k-1] = order[k,i];
                        
                    }
                    Vector mparmtemp = (Vector)mparm[":", i];
                    Vector tparmtemp = (Vector)tparm[":", i];
                    parity(orderKarray, nb[1,i], ref mparmtemp, ref tparmtemp);
                    tparm[":", i] = Matrix<double>.Transpose(tparmtemp);
                    mparm[":", i] = Matrix<double>.Transpose(mparmtemp);
                }
                // Set transition levels

                alttns = vnode[nlev - 1];
                altsym = vnode[nlev];
                altiso = vnode[nlev+1];

                // Allocate the global store of quasi-static parameters

                maxo = Math.Max(maxs, Math.Max(maxm, maxl));
                omaxhwm = maxo;
                nmaxhwm = maxn;

                fs = new Matrix<double>(maxs + 1, 2);
                fl = new Matrix<double>(maxl + 1, 2);
                fm = new Matrix<double>(maxm + 1, 2);
                bz = new Vector(nbf);
                bm = new Vector(nbf);
                zwght = new Vector(p + 1);

                // change the initalization flag and reset some other things

                previous = new Vector(new List<double>(new double[] { -1, -1, -1, -1, -1 }));
                //parent = hwm;

            }
            public Vector hwmqt(int IYD, double SEC, double ALT, double GLAT, double GLON, double STL, double F107A, double F107)
            {
                // ------------------------------------------------------------
                // The quiet time only HWM function call
                // ------------------------------------------------------------

                // Local variables
                double deg2rad =  Math.PI/180;
                Vector input = new Vector(5);
                double u, v;
                double cs, ss, cm, sm, cl, sl;
                double AA, BB, CC, DD;
                double vb, wb;
                double theta, sc;

                int b, c, d, m, n, s, l;

                int amaxs, amaxn;
                int pmaxm, pmaxs, pmaxn;
                int tmaxl, tmaxs, tmaxn;

                bool[] refresh = new bool[5];

                // ====================================================================
                // Update VSH model terms based on any change in the input parameters
                // ====================================================================

                //if (qwminit)
                //initqwm(qwmdefault);

                //IYD = IYD % 1000;
                input[1] = IYD % 1000;
                input[2] = SEC;
                input[3] = GLON;
                input[4] = GLAT;
                input[5] = ALT;

                // Seasonal variations
                if (input[1] != previous[1])
                {
                    AA = input[1] * 2 * Math.PI / 365.25;
                    for (s = 1; s <= maxs + 1; s++)
                    {
                        BB = (s - 1) * AA;
                        fs[s, 1] = Math.Cos(BB);
                        fs[s, 2] = Math.Sin(BB);
                    }
                    refresh = Array.ConvertAll(refresh, x => true);
                    previous[1] = input[1];
                }

                // Hourly time changes, tidal variations

                if (input[2] != previous[2] || input[3] != previous[3])
                {
                    AA = (input[2] / 3600 + input[3] / 15 + 48) % 24;
                    BB = AA * 2 * Math.PI / 24;
                    for (l = 1; l <= maxl + 1; l++)
                    {
                        CC = (l - 1) * BB;
                        fl[l, 1] = Math.Cos(CC);
                        fl[l, 2] = Math.Sin(CC);
                    }
                    refresh[3] = true;   // tides
                    previous[2] = input[2];
                }

                // Longitudinal variations, stationary planetary waves

                if (input[3] != previous[3])
                {
                    AA = input[3] * deg2rad;
                    for (m = 1; m <= maxm + 1; m++)
                    {
                        BB = (m - 1) * AA;
                        fm[m, 1] = Math.Cos(BB);
                        fm[m, 2] = Math.Sin(BB);
                    }
                    refresh[2] = true;   // stationary planetary waves
                    previous[3] = input[3];
                }

                // Latitude

                theta = (90.0 - input[4]) * deg2rad;
                if (input[4] != glatalf)
                {
                    AA = (90.0 - input[4]) * deg2rad;    // theta = colatitude in radians
                    alfbasis(maxn, maxm, AA, out gpbar, out gvbar, out gwbar);
                    refresh[1] = true;
                    refresh[2] = true;
                    refresh[3] = true;
                    refresh[4] = true;
                    glatalf = input[4];
                    previous[4] = input[4];
                }

                // Altitude

                if (input[5] != previous[5])
                {
                    vertwght(input[5], out zwght, out lev);
                    previous[5] = input[5];
                }

                // ====================================================================
                // Calculate the VSH functions
                // ====================================================================

                u = 0.0;
                v = 0.0;

                for (b = 1; b <= p+1; b++)
                {
                    if (zwght[b] == 0) 
                        continue;

                    d = b + lev-1;

                    if (priornb != nb[1,d])
                        refresh.All(x => true); // recalculate basis functions
                    priornb = nb[1,d];

                    if (!refresh.All(x => x))
                    {
                        c = nb[1,d];
                        if (component[0])
                                u = u + zwght[b] * Vector.Dot(bz[new MatrixIndex(1,c)], (Vector)mparm[new MatrixIndex(1, c), d]);
                        if (component[1])
                                v = v + zwght[b] * Vector.Dot(bz[new MatrixIndex(1, c)], (Vector)tparm[new MatrixIndex(1, c), d]);
                        continue;
                    }

                    amaxs = order[1, d];
                    amaxn = order[2, d];
                    pmaxm = order[3, d];
                    pmaxs = order[4, d];
                    pmaxn = order[5, d];
                    tmaxl = order[6, d];
                    tmaxs = order[7, d];
                    tmaxn = order[8, d];

                    c = 1;

                    // ------------- Seasonal - Zonal average(m = 0) ----------------

                    if (refresh[1] && content[1])
                    {
                        for (n = 1; n <= amaxn; n++)               // s = 0
                        {
                            bz[c] = -Math.Sin(n * theta);   //
                            bz[c + 1] = Math.Sin(n * theta);
                            c = c + 2;
                        }
                        for (s = 1; s <= amaxs; s++) // Seasonal variations
                        {
                            cs = fs[s+1, 1];
                            ss = fs[s+1, 2];
                            for (n = 1; n <= amaxn; n++)
                            {
                                sc = Math.Sin(n * theta);
                                bz[c] = -sc * cs;   // Cr A
                                bz[c + 1] = sc * ss;  // Ci B
                                bz[c + 2] = sc * cs;
                                bz[c + 3] = -sc * ss;
                                c = c + 4;
                            }
                        }
                        cseason = c;
                    }
                    else
                        c = cseason;

                    // ---------------- Stationary planetary waves --------------------

                    if (refresh[2] && content[2])
                    {
                        for (m = 1; m <= pmaxm; m++)
                        {
                            cm = fm[m+1, 1] * wavefactor[m+1];
                            sm = fm[m+1, 2] * wavefactor[m+1];
                            for (n = m; n <= pmaxn; n++)
                            {// s = 0
                                vb = gvbar[n+1, m+1];
                                wb = gwbar[n+1, m+1];
                                bz[c] = -vb * cm;    // Cr* (cm) * -vb A
                                bz[c + 1] = vb * sm;    // Ci* (sm) * vb   B
                                bz[c + 2] = -wb * sm;   // Br* (sm) * -wb C
                                bz[c + 3] = -wb * cm;   // Bi* (cm) * -wb D
                                c = c + 4;
                            }
                            for (s = 1; s <= pmaxs; s++) 
                            {
                                cs = fs[s+1, 1];
                                ss = fs[s+1, 2];
                                for (n = m; n <= pmaxn; n++)
                                {
                                    vb = gvbar[n+1, m+1];
                                    wb = gwbar[n+1, m+1];
                                    bz[c] = -vb * cm * cs;  // Crc* (cmcs) * -vb A
                                    bz[c + 1] = vb * sm * cs; // Cic* (smcs) * vb   B
                                    bz[c + 2] = -wb * sm * cs;  // Brc* (smcs) * -wb C
                                    bz[c + 3] = -wb * cm * cs;  // Bic* (cmcs) * -wb D
                                    bz[c + 4] = -vb * cm * ss;  // Crs* (cmss) * -vb E
                                    bz[c + 5] = vb * sm * ss; // Cis* (smss) * vb   F
                                    bz[c + 6] = -wb * sm * ss;  // Brs* (smss) * -wb G
                                    bz[c + 7] = -wb * cm * ss;  // Bis* (cmss) * -wb H
                                    c = c + 8;
                                }

                            }
                            cwave = c;
                        }
                    }
                    else
                        c = cwave;

                    // ---------------- Migrating Solar Tides ---------------------

                    if (refresh[3] && content[3])
                    {
                        for (l = 1; l <= tmaxl; l++)
                        {
                            cl = fl[l+1, 1] * tidefactor[l+1];
                            sl = fl[l+1, 2] * tidefactor[l+1];
                            for (n = l; n <= tmaxn; n++)           // s = 0
                            {
                                vb = gvbar[n+1, l+1];
                                wb = gwbar[n+1, l+1];
                                bz[c] = -vb * cl;    // Cr* (cl) * -vb
                                bz[c + 1] = vb * sl;   // Ci* (sl) * vb
                                bz[c + 2] = -wb * sl;   // Br* (sl) * -wb
                                bz[c + 3] = -wb * cl;    // Bi* (cl) * -wb
                                c = c + 4;
                            }
                            for (s = 1; s <= tmaxs; s++)
                            {
                                cs = fs[s+1, 1];
                                ss = fs[s+1, 2];
                                for (n = l; n <= tmaxn; n++)
                                {
                                    vb = gvbar[n+1, l+1];
                                    wb = gwbar[n+1, l+1];
                                    bz[c] = -vb * cl * cs;  // Crc* (clcs) * -vb
                                    bz[c + 1] = vb * sl * cs; // Cic* (slcs) * vb
                                    bz[c + 2] = -wb * sl * cs;  // Brc* (slcs) * -wb
                                    bz[c + 3] = -wb * cl * cs;   // Bic* (clcs) * -wb
                                    bz[c + 4] = -vb * cl * ss;   // Crs* (clss) * -vb
                                    bz[c + 5] = vb * sl * ss; // Cis* (slss) * vb
                                    bz[c + 6] = -wb * sl * ss;   // Brs* (slss) * -wb
                                    bz[c + 7] = -wb * cl * ss;   // Bis* (clss) * -wb
                                    c = c + 8;
                                }
                            }
                            ctide = c;
                        }
                    }
                    else
                        c = ctide;

                    // ---------------- Non-Migrating Solar Tides ------------------

                    // TBD

                    c = c - 1;

                    // ====================================================================
                    // Calculate the wind components
                    // ====================================================================

                    if (component[0])
                    {
                            u = u + zwght[b] * Matrix<double>.Dot(bz[new MatrixIndex(1, c)], mparm[new MatrixIndex(1, c), d]);
                    }
                    if (component[1])
                    {
                            v = v + zwght[b] * Matrix<double>.Dot(bz[new MatrixIndex(1,c)], tparm[new MatrixIndex(1, c), d]);
                    }

                }
                Vector w = new Vector(2);
                w[1] = v;
                w[2] = u;
                return w;

            }
            private void vertwght(double alt, out Vector wght, out int iz)
            {

                //real[8],intent(in)      :: alt
                //real[8],intent(out)     :: wght[4]
                //integer[4],intent(out)  :: iz

                //real[8]             :: we(0:4)
                Vector we = new Vector(5);
                iz = 0;
                wght = new Vector(4);
                iz = findspan(nnode - p - 1, p, alt, vnode) - p;
                iz = Math.Min(iz, 27);

                wght[1] = bspline(p, nnode, vnode, iz, alt);
                wght[2] = bspline(p, nnode, vnode, iz + 1, alt);
                if (iz <= 26)
                {
                    wght[3] = bspline(p, nnode, vnode, iz + 2, alt);
                    wght[4] = bspline(p, nnode, vnode, iz + 3, alt);
                    return;
                }
                if (alt > alttns)
                {
                    we[1] = 0.0;
                    we[2] = 0.0;
                    we[3] = 0.0;
                    we[4] = Math.Exp(-(alt - alttns) / H);
                    we[5] = 1.0;
                }
                else
                {
                    we[1] = bspline(p, nnode, vnode, iz + 2, alt);
                    we[2] = bspline(p, nnode, vnode, iz + 3, alt);
                    we[3] = bspline(p, nnode, vnode, iz + 4, alt);
                    we[4] = 0.0;
                    we[5] = 0.0;
                }
                wght[3] = Matrix<double>.Dot(we, e1);
                wght[4] = Matrix<double>.Dot(we, e2);

                return;
            }
            public static double bspline(int p, int m, Vector V, int i, double u)
            {
                Vector N = new Vector(p + 2);
                double Vleft, Vright;
                double saved, temp;
                int j, k;

                if ((i == 0) && (u == V[1]))
                    return 1.0;
                if ((i == (m - p - 1)) && (u == V[m+1]))
                    return 1.0;

                if (u < V[i] || u >= V[i + p + 1])
                    return 0.0;

                //N = 0.0;
                for (j = 1; j <= p+1; j++)
                {
                    if (u >= V[i + j - 1] && u < V[i + j])
                        N[j] = 1.0;
                    else
                        N[j] = 0.0;
                }

                for (k = 1; k <= p; k++)
                {
                    if (N[1] == 0)
                        saved = 0;
                    else
                        saved = ((u - V[i]) * N[1]) / (V[i + k] - V[i]);
                    for (j = 1; j <= p - k+1; j++)
                    {
                        Vleft = V[i + j];
                        Vright = V[i + j + k];
                        if (N[j+1] == 0)
                        {
                            N[j] = saved;
                            saved = 0;
                        }
                        else
                        {
                            temp = N[j+1] / (Vright - Vleft);
                            N[j] = saved + (Vright - u) * temp;
                            saved = (u - Vleft) * temp;
                        }
                    }
                }
                return N[1];

            }
            private int findspan(int n, int p, double u, Vector V)
            {
                // =====================================================
                // Function to locate the knot span
                // =====================================================
                int low, mid, high;

                if (u >= V[n + 1])
                    return n;

                low = p;
                high = n + 1;
                mid = (low + high) / 2+1;

                while (u < V[mid] || u >= V[mid + 1])
                {
                    if (u < V[mid])
                        high = mid;
                    else
                        low = mid;
                    mid = (low + high) / 2;
                }


                return mid;
            }
            private void parity(int[] order, int nb, ref Vector mparm, ref Vector tparm)
            {
                //integer[4],intent(in)     :: order[8]
                //integer[4],intent(in)     :: nb
                //real[8],intent(inout)     :: mparm(nb)
                //real[8],intent(out)       :: tparm(nb)

                int c, m, n, s, l;

                int amaxs, amaxn;
                int pmaxm, pmaxs, pmaxn;
                int tmaxl, tmaxs, tmaxn;
                //tparm = new Vector(nb);
                //mparm = new Vector(nb);
                amaxs = order[0];
                amaxn = order[1];
                pmaxm = order[2];
                pmaxs = order[3];
                pmaxn = order[4];
                tmaxl = order[5];
                tmaxs = order[6];
                tmaxn = order[7];

                c = 1;

                for (n = 1; n <= amaxn; n++)
                {
                    tparm[c] = 0.0;
                    tparm[c + 1] = -mparm[c + 1];
                    mparm[c + 1] = 0.0;
                    c = c + 2;
                }
                for (s = 1; s <= amaxs; s++)
                {
                    for (n = 1; n <= amaxn; n++)
                    {
                        tparm[c] = 0.0;
                        tparm[c + 1] = 0.0;
                        tparm[c + 2] = -mparm[c + 2];
                        tparm[c + 3] = -mparm[c + 3];
                        mparm[c + 2] = 0.0;
                        mparm[c + 3] = 0.0;
                        c = c + 4;
                    }
                }

                for (m = 1; m <= pmaxm; m++)
                {
                    for (n = m; n <= pmaxn; n++)
                    {
                        tparm[c] = mparm[c + 2];
                        tparm[c + 1] = mparm[c + 3];
                        tparm[c + 2] = -mparm[c];
                        tparm[c + 3] = -mparm[c + 1];
                        c = c + 4;
                    }
                    for (s = 1; s <= pmaxs; s++)
                    {
                        for (n = m; n <= pmaxn; n++)
                        {
                            tparm[c] = mparm[c + 2];
                            tparm[c + 1] = mparm[c + 3];
                            tparm[c + 2] = -mparm[c];
                            tparm[c + 3] = -mparm[c + 1];
                            tparm[c + 4] = mparm[c + 6];
                            tparm[c + 5] = mparm[c + 7];
                            tparm[c + 6] = -mparm[c + 4];
                            tparm[c + 7] = -mparm[c + 5];
                            c = c + 8;
                        }
                    }

                }
                for (l = 1; l <= tmaxl; l++)
                {
                    for (n = l; n <= tmaxn; n++)
                    {
                        tparm[c] = mparm[c + 2];
                        tparm[c + 1] = mparm[c + 3];
                        tparm[c + 2] = -mparm[c];
                        tparm[c + 3] = -mparm[c + 1];
                        c = c + 4;
                    }
                    for (s = 1; s <= tmaxs; s++)
                    {
                        for (n = l; n <= tmaxn; n++)
                        {
                            tparm[c] = mparm[c + 2];
                            tparm[c + 1] = mparm[c + 3];
                            tparm[c + 2] = -mparm[c];
                            tparm[c + 3] = -mparm[c + 1];
                            tparm[c + 4] = mparm[c + 6];
                            tparm[c + 5] = mparm[c + 7];
                            tparm[c + 6] = -mparm[c + 4];
                            tparm[c + 7] = -mparm[c + 5];
                            c = c + 8;
                        }
                    }
                }

                //return tparm;

            }
        }
        public class DWM
        {

            int nterm;             // Number of terms in the model
            int nmax, mmax;         // Max latitudinal degree
            int nvshterm;         // # of VSH basis functions

            Matrix<int> termarr;      // 3 x nterm index of coupled terms
            Vector coeff;         // Model coefficients
            Matrix<double> vshterms;     // VSH basis values
            Matrix<double> termval;      // Term values to which coefficients are applied
            Matrix<double> dpbar;        // Associated lengendre fns
            Matrix<double> dvbar;
            Matrix<double> dwbar;
            Matrix<double> mltterms;    // MLT Fourier terms
            double twidth;            // Transition width of high-lat mask

            public static double sineps = 0.39781868;
            double deg2rad = Math.PI / 180;
            double mltlast = 1e16, mlatlast = 1e16, kplast = 1e16;
            double glatlast = 1.0e16, glonlast = 1.0e16;
            double daylast = 1.0e16, utlast = 1.0e16, aplast = 1.0e16;
            double day, ut, mlat, mlon, mlt, kp;
            double f1e, f1n, f2e, f2n;
            double talt = 125.0;
            Vector kpterms = new Vector(3);
            gd2qdc gd;
            //HWM parent;

            public DWM(out int nmaxout, out int mmaxout)
            {
                nterm = HWMData.dwm.nterm;
                mmax = HWMData.dwm.mmax;
                nmax = HWMData.dwm.nmax;
                termarr = HWMData.dwm.termarr;
                coeff = HWMData.dwm.coeff;
                twidth = HWMData.dwm.twidth;
                nvshterm = (((nmax + 1) * (nmax + 2) - (nmax - mmax) * (nmax - mmax + 1)) / 2 - 1) * 4 - 2 * nmax;
                termval = new Matrix<double>(2, nterm);
                dpbar = new Matrix<double>(nmax + 1, mmax + 1);
                dvbar = new Matrix<double>(nmax + 1, mmax + 1);
                dwbar = new Matrix<double>(nmax + 1, mmax + 1);
                mltterms = new Matrix<double>(mmax + 1, 2);
                vshterms = new Matrix<double>(2, nvshterm);
                nmaxout = nmax;
                mmaxout = mmax;


                gd = new gd2qdc();

                //parent = hwm;
            }
            public Vector dwm07(int iyd, double sec, double alt, double glat, double glon, Vector ap)
            {
                Vector dw = new Vector(2);
                
                double mmpwind, mzpwind;
                
                
                

                //CONVERT AP TO KP
                if (ap[2] != aplast)
                    kp = ap2kp(ap[2]);


                //CONVERT GEO LAT/LON TO QD LAT/LON
                if ((glat != glatlast) || (glon != glonlast))
                    gd.gd2qd(glat, glon, out mlat, out mlon, out f1e, out f1n, out f2e, out f2n);


                //COMPUTE QD MAGNETIC LOCAL TIME(LOW-PRECISION)
                day = iyd % 1000;
                ut = sec / 3600.0;
                if ((day != daylast) || (ut != utlast) || (glat != glatlast) || (glon != glonlast))
                    mlt = gd.mltcalc(mlat, mlon, day, ut);

                //RETRIEVE DWM WINDS
                Dwm07b(mlt, mlat, kp, out mmpwind, out mzpwind);

                //CONVERT TO GEOGRAPHIC COORDINATES
                dw[1] = f2n * mmpwind + f1n * mzpwind;
                dw[2] = f2e * mmpwind + f1e * mzpwind;

                //APPLY HEIGHT PROFILE
                dw = dw / (1 + Math.Exp(-(alt - talt) / twidth));

                glatlast = glat;
                glonlast = glon;
                daylast = day;
                utlast = ut;
                aplast = ap[2];

                return dw;
            }
            private void Dwm07b(double mlt, double mlat, double kp, out double mmpwind, out double mzpwind)
            {

                // Local variables
                int iterm, ivshterm, n, m;
                Vector termvaltemp;
                double latwgtterm;
                double theta, phi, mphi;
                mmpwind = 0;
                mzpwind = 0;

                //COMPUTE LATITUDE PART OF VSH TERMS
                if (mlat != mlatlast)
                {
                    theta = (90 - (mlat)) * deg2rad;
                    alfbasis(nmax, mmax, theta, out dpbar, out dvbar, out dwbar);
                }

                //COMPUTE MLT PART OF VSH TERMS
                if (mlt != mltlast)
                {
                    phi = (mlt) * deg2rad * 15;
                    for (m = 1; m <= mmax + 1; m++)
                    {
                        mphi = (m - 1) * phi;
                        mltterms[m, 1] = Math.Cos(mphi);
                        mltterms[m, 2] = Math.Sin(mphi);
                    }
                }

                //COMPUTE VSH TERMS
                if ((mlat != mlatlast) || (mlt != mltlast))
                {
                    ivshterm = 1;
                    for (n = 2; n <= nmax+1; n++)
                    {
                        vshterms[1, ivshterm] = -(dvbar[n, 1] * mltterms[1, 1]);
                        vshterms[1, ivshterm + 1] = (dwbar[n, 1] * mltterms[1, 1]);
                        vshterms[2, ivshterm] = -vshterms[1, ivshterm + 1];
                        vshterms[2, ivshterm + 1] = vshterms[1, ivshterm];
                        ivshterm = ivshterm + 2;
                        for (m = 2; m <= mmax+1; m++)
                        {
                            if (m > n)
                                continue;
                            vshterms[1, ivshterm] = -(dvbar[n, m] * mltterms[m, 1]);
                            vshterms[1, ivshterm + 1] = (dvbar[n, m] * mltterms[m, 2]);
                            vshterms[1, ivshterm + 2] = (dwbar[n, m] * mltterms[m, 2]);
                            vshterms[1, ivshterm + 3] = (dwbar[n, m] * mltterms[m, 1]);
                            vshterms[2, ivshterm] = -vshterms[1, ivshterm + 2];
                            vshterms[2, ivshterm + 1] = -vshterms[1, ivshterm + 3];
                            vshterms[2, ivshterm + 2] = vshterms[1, ivshterm];
                            vshterms[2, ivshterm + 3] = vshterms[1, ivshterm + 1];
                            ivshterm = ivshterm + 4;
                        }
                    }
                }

                //COMPUTE KP TERMS
                if (kp != kplast)
                    kpterms = kpspl3(kp);
                    

                //COMPUTE LATITUDINAL WEIGHTING TERM
                latwgtterm = latwgt2(mlat, mlt, kp, twidth);

                //GENERATE COUPLED TERMS
                for (iterm = 1; iterm <= nterm; iterm++)
                {
                    termvaltemp = new Vector(new List<double>(new double[] { 1.0, 1.0 }));
                    if (termarr[1, iterm] != 999)
                    {
                        termvaltemp[1] = termvaltemp[1] * vshterms[1, termarr[1, iterm]+1];
                        termvaltemp[2] = termvaltemp[2] * vshterms[2, termarr[1, iterm]+1];
                    }
                    if (termarr[2, iterm] != 999)
                        termvaltemp = termvaltemp * kpterms[termarr[2, iterm]+1];
                    if (termarr[3, iterm] != 999)
                        termvaltemp = termvaltemp * latwgtterm;
                    termval[1, iterm] = termvaltemp[1];
                    termval[2, iterm] = termvaltemp[2];
                }


                //APPLY COEFFICIENTS
                mmpwind = Matrix<double>.Dot(coeff, termval[1, new MatrixIndex(1,nterm)]);
                mzpwind = Matrix<double>.Dot(coeff, termval[2, new MatrixIndex(1, nterm)]);
                
                mlatlast = mlat;
                mltlast = mlt;
                kplast = kp;
            }
            
            private Vector kpspl3(double kp)
            {
                //================================================================================
                //                           Cubic Spline interpolation of Kp
                //===============================================================================
                Vector kpterms = new Vector(3);
                int i, j;
                double x;
                double[] kpspl = new double[7];
                double[] node = new double[] { -10, -8, 0, 2, 5, 8, 18, 20 };

                x = Math.Max(kp, 0.0);
                x = Math.Min(x, 8.0);

                for (i = 0; i <= 6; i++)
                {
                    kpspl[i] = 0.0;
                    if ((x >= node[i]) && (x < node[i + 1]))
                        kpspl[i] = 1.0;
                }
                for (j = 2; j <= 3; j++)
                {
                    for (i = 0; i <= 8 - j - 1; i++)
                    {
                        kpspl[i] = kpspl[i] * (x - node[i]) / (node[i + j - 1] - node[i]) + kpspl[i + 1] * (node[i + j] - x) / (node[i + j] - node[i + 1]);
                    }
                }
                kpterms[1] = kpspl[0] + kpspl[1];
                kpterms[2] = kpspl[2];
                kpterms[3] = kpspl[3] + kpspl[4];

                return kpterms;

            }
            private double latwgt2(double mlat, double mlt, double kp0, double twidth)
            {   //================================================================================
                //                           (Function) Latitude weighting factors
                //================================================================================

                double kp;
                double mltrad, sinmlt, cosmlt, tlat;

                double[] coeff = new double[] { 65.7633, -4.60256, -3.53915, -1.99971, -0.752193, 0.972388 };


                mltrad = mlt * 15.0 * deg2rad;
                sinmlt = Math.Sin(mltrad);
                cosmlt = Math.Cos(mltrad);
                kp = Math.Max(kp0, 0.0);
                kp = Math.Min(kp, 8.0);
                tlat = coeff[0] + coeff[1] * cosmlt + coeff[2] * sinmlt + kp * (coeff[3] + coeff[4] * cosmlt + coeff[5] * sinmlt);
                return 1.0 / (1 + Math.Exp(-(Math.Abs(mlat) - tlat) / twidth));
            }
            private double ap2kp(double ap0)
            {    //=================================================================================
                 //                           Convert Ap to Kp
                 //=================================================================================

                Vector apgrid = new Vector(new List<double>(new double[] {0 ,2 ,3 ,4 ,5 ,6 ,7 ,9 ,12 ,15 ,18 ,
                                                                     22 ,27 ,32 ,39 ,48 ,56 ,67 ,80 ,94 ,
                                                                     111 ,132 ,154 ,179 ,207 ,236 ,300 ,400 }));
                Vector kpgrid = new Vector(new List<double>(new double[]{0 ,1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,
                                                                     12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,
                                                                     22 ,23 ,24 ,25 ,26 ,27 })) / 3;
                double ap;
                int i;


                ap = ap0;
                if (ap < 0)
                    ap = 0;
                if (ap > 400)
                    ap = 400;

                i = 1;
                while (ap > apgrid[i])
                    i = i + 1;
                if (ap == apgrid[i])
                    return kpgrid[i];
                else
                    return kpgrid[i - 1] + (ap - apgrid[i - 1]) / (3.0 * (apgrid[i] - apgrid[i - 1]));
            }
        }
        public class gd2qdc
        {
            // ########################################################################
            //     Geographic <=> Geomagnetic Coordinate Transformations
            //
            //  Converts geodetic coordinates to Quasi-Dipole coordinates(Richmond, J.Geomag.
            //  Geoelec., 1995, p. 191), using a spherical harmonic representation.
            //
            // ########################################################################

            static int nterm, nmax, mmax;  //Spherical harmonic expansion parameters

            Matrix<double> coeff;      //Coefficients for spherical harmonic expansion

            static Matrix<double> shgradtheta;    //Array to hold spherical harmonic gradients
            static Matrix<double> shgradphi;    //Array to hold spherical harmonic gradients
            static Matrix<double> normadj;   //Adjustment to VSH normalization factor
            static double epoch, alt;

            static double deg2rad = Math.PI / 180.0;

            //HWM parent;

            public gd2qdc()
            {

                int iterm, n;
                int j;

                nmax = HWMData.gd2qd.nmax;
                mmax = HWMData.gd2qd.mmax;
                nterm = HWMData.gd2qd.nterm;
                epoch = HWMData.gd2qd.epoch;
                alt = HWMData.gd2qd.alt;
                coeff = HWMData.gd2qd.coeff;



                xcoeff = new Matrix<double>(1, nterm);
                ycoeff = new Matrix<double>(1, nterm);
                zcoeff = new Matrix<double>(1, nterm);
                sh = new Matrix<double>(1, nterm);
                shgradtheta = new Matrix<double>(1, nterm);
                shgradphi = new Matrix<double>(1, nterm);
                normadj = new Matrix<double>(1, nmax + 1);


                for (iterm = 1; iterm <= nterm; iterm++)
                {
                    xcoeff[iterm] = coeff[iterm, 1];
                    ycoeff[iterm] = coeff[iterm, 2];
                    zcoeff[iterm] = coeff[iterm, 3];
                }

                for (n = 0; n <= nmax; n++)
                    normadj[n+1] = Math.Sqrt((n * (n + 1)));


                nmaxqdc = nmax;
                mmaxqdc = mmax;

                //parent = hwm;
            }

            public void gd2qd(double glat, double glon, out double qlat, out double qlon, out double f1e, out double f1n, out double f2e, out double f2n)
            {
                int n, m, i;
                double theta, phi;
                double mphi, cosmphi, sinmphi;
                double x, y, z;
                double cosqlat, cosqlon, sinqlon;
                double xgradtheta, ygradtheta, zgradtheta;
                double xgradphi, ygradphi, zgradphi;
                double qlonrad;



                if (glat != glatalf)
                {
                    theta = (90 - glat) * deg2rad;
                    alfbasis(nmax, mmax, theta, out gpbar, out gvbar, out gwbar);
                    glatalf = glat;
                }
                phi = glon * deg2rad;


                i = 1;
                for (n = 0; n <= nmax; n++)
                {
                    sh[i] = gpbar[n+1, 1];
                    shgradtheta[i] = gvbar[n+1, 1] * normadj[n+1];
                    shgradphi[i] = 0;
                    i = i + 1;
                }
                for (m = 1; m <= mmax; m++)
                {
                    mphi = m * phi;
                    cosmphi = Math.Cos(mphi);
                    sinmphi = Math.Sin(mphi);
                    for (n = m; n <= nmax; n++)
                    {
                        sh[i] = gpbar[n+1, m + 1] * cosmphi;
                        sh[i + 1] = gpbar[n + 1, m + 1] * sinmphi;
                        shgradtheta[i] = gvbar[n + 1, m + 1] * normadj[n + 1] * cosmphi;
                        shgradtheta[i + 1] = gvbar[n + 1, m + 1] * normadj[n + 1] * sinmphi;
                        shgradphi[i] = -gwbar[n + 1, m + 1] * normadj[n + 1] * sinmphi;
                        shgradphi[i + 1] = gwbar[n + 1, m + 1] * normadj[n + 1] * cosmphi;
                        i = i + 2;
                    }
                }



                x = Matrix<double>.Dot(sh, xcoeff);
                y = Matrix<double>.Dot(sh, ycoeff);
                z = Matrix<double>.Dot(sh, zcoeff);


                qlonrad = Math.Atan2(y, x);
                cosqlon = Math.Cos(qlonrad);
                sinqlon = Math.Sin(qlonrad);
                cosqlat = x * cosqlon + y * sinqlon;


                qlat = Math.Atan2(z, cosqlat) / deg2rad;
                qlon = qlonrad / deg2rad;


                xgradtheta = Matrix<double>.Dot(shgradtheta, xcoeff);
                ygradtheta = Matrix<double>.Dot(shgradtheta, ycoeff);
                zgradtheta = Matrix<double>.Dot(shgradtheta, zcoeff);


                xgradphi = Matrix<double>.Dot(shgradphi, xcoeff);
                ygradphi = Matrix<double>.Dot(shgradphi, ycoeff);
                zgradphi = Matrix<double>.Dot(shgradphi, zcoeff);


                f1e = -zgradtheta * cosqlat + (xgradtheta * cosqlon + ygradtheta * sinqlon) * z;
                f1n = -zgradphi * cosqlat + (xgradphi * cosqlon + ygradphi * sinqlon) * z;
                f2e = ygradtheta * cosqlon - xgradtheta * sinqlon;
                f2n = ygradphi * cosqlon - xgradphi * sinqlon;

            }

            public double mltcalc(double qlat, double qlon, double day, double ut)
            {
                //==================================================================================
                //                  (Function) Calculate Magnetic Local Time
                //==================================================================================
                int n, m, i;
                double asunglat, asunglon, asunqlon;
                double glat, theta, phi;
                double mphi, cosmphi, sinmphi;
                double x, y;
                double cosqlat, cosqlon, sinqlon;
                double qlonrad;

                //COMPUTE GEOGRAPHIC COORDINATES OF ANTI-SUNWARD DIRECTION(LOW PRECISION)
                asunglat = -Math.Asin(Math.Sin((day + ut / 24.0 - 80.0) * deg2rad) * DWM.sineps) / deg2rad;
                asunglon = -ut * 15;

                //COMPUTE MAGNETIC COORDINATES OF ANTI-SUNWARD DIRECTION
                theta = (90 - asunglat) * deg2rad;
                alfbasis(nmax, mmax, theta, out spbar, out svbar, out swbar);
                phi = asunglon * deg2rad;
                i = 1;
                for (n = 0; n <= nmax; n++)
                {
                    sh[i] = spbar[n + 1, 1];
                    i = i + 1;
                }
                for (m = 1; m <= mmax; m++)
                {
                    mphi = m * phi;
                    cosmphi = Math.Cos(mphi);
                    sinmphi = Math.Sin(mphi);
                    for (n = m; n <= nmax; n++)
                    {
                        sh[i] = spbar[n + 1, m + 1] * cosmphi;
                        sh[i + 1] = spbar[n + 1, m + 1] * sinmphi;
                        i = i + 2;
                    }
                }
                x = Matrix<double>.Dot(sh, xcoeff);
                y = Matrix<double>.Dot(sh, ycoeff);
                asunqlon = Math.Atan2(y, x) / deg2rad;

                //COMPUTE MLT
                return (qlon - asunqlon) / 15.0;
            }
        }
    }
}