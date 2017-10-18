using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Collections;

using GH_IO;
using GH_IO.Serialization;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using System;
using System.IO;
using System.Xml;
using System.Xml.Linq;
using System.Linq;
using System.Data;
using System.Drawing;
using System.Reflection;
using System.Collections;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using IronPython;

using Microsoft;

namespace BIG
{
    public class RecolorMesh : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public RecolorMesh()
          : base("RecolorMesh", "RecolorMesh",
              "Recoloring mesh with parallel threading",
              "BIG", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "Mesh", "Mesh to recolor", GH_ParamAccess.item);
            pManager.AddNumberParameter("Result", "Result", "Result with the same number of faces as mesh", GH_ParamAccess.list);
            pManager.AddColourParameter("List Of Colors", "Colors", "Optional list of colours as gradients", GH_ParamAccess.list);
            pManager.AddNumberParameter("Lower Boundary", "LowBound", "Optional lower bound for the coloring", GH_ParamAccess.item);
            pManager.AddNumberParameter("Upper Boundary", "UpperBound", "Optional upper bound for the coloring", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Recolored", "NewMesh", "Recoloured mesh based on result and colors", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //1.0 Collecting data

            Mesh mesh = new Mesh();
            List<double> result = new List<double>();
            List<Color> coloraslist = new List<Color>();

            double ming = 0.0;
            double maxg = 0.0;

            //1.1 Return conditions
            if ((!DA.GetData(0, ref mesh)))
                return;

            if ((!DA.GetDataList(1, result)))
                return;

            if (!DA.GetDataList(2, coloraslist))
                return;

            if (!DA.GetData(3, ref ming))
                return;

            if (!DA.GetData(4, ref maxg))
                return;


            //2.0 Setting up the run;

            if (mesh.Faces.Count() != result.Count()) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "unequal faces and result"); return; }

            List<Mesh> meshall = new List<Mesh>();
            Mesh ms = new Mesh();


            var numbers = Enumerable.Range(0, mesh.Faces.Count());
            var _result = numbers.AsParallel().AsOrdered();

            int fA = 0;

            foreach (var i in _result)
            {
                Color cf = gradientlist(coloraslist.ToArray(), ming, maxg).ColourAt(result[i]);
                MeshFace face = mesh.Faces[i];

                ms.Vertices.Add(mesh.Vertices[face.A]);
                ms.Vertices.Add(mesh.Vertices[face.B]);
                ms.Vertices.Add(mesh.Vertices[face.C]);

                ms.VertexColors.Add(cf);
                ms.VertexColors.Add(cf);
                ms.VertexColors.Add(cf);

                if (face.IsQuad)
                {
                    ms.Vertices.Add(mesh.Vertices[face.D]);
                    ms.VertexColors.Add(cf);
                    ms.Faces.AddFace(fA, fA + 1, fA + 2, fA + +3);
                    fA = fA + 4;

                }
                else
                {
                    ms.Faces.AddFace(fA, fA + 1, fA + 2);
                    fA = fA + 3;
                }
            }
            DA.SetData(0, ms);
        }

        // X. Extra additional useful functions
        // X.01 Gradient maker 
        private Grasshopper.GUI.Gradient.GH_Gradient gradientlist(Color[] colorarray, double t0, double t1)
        {
            Grasshopper.GUI.Gradient.GH_Gradient gradient2 = new Grasshopper.GUI.Gradient.GH_Gradient();
            for (int i = 0; i < colorarray.Count(); i++)
            {
                double grip = t0 + (double)i * (t1 - t0) / (colorarray.Count() - 1); //fix 10.10.17 add t0, to rescale the gradient grip
                gradient2.AddGrip(grip, colorarray[i]);
            }
            return gradient2;
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                return Properties.Resources.recolormesh;
                //return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("0cdcf43d-d82e-42e3-a55f-79cb4a67cafd"); }
        }
    }

    public class SolarPath : GH_Component
    {
        public SolarPath()
            : base("Solar Path", "SolarPath",
              "Produce only vectors of a given Location on Earth and timesteps",
              "BIG", "Analysis")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddVectorParameter("North", "North", "Direction of north", GH_ParamAccess.item, Vector3d.YAxis); //DA.0
            pManager.AddNumberParameter("Timezone", "TimeZone", "Timezone of the area", GH_ParamAccess.item);//DA.1
            pManager.AddNumberParameter("Latitude", "Lat", "Latitude of the location coordinates", GH_ParamAccess.item);//DA.2
            pManager.AddNumberParameter("Longitude", "Lon", "Longitude of the location coordinates", GH_ParamAccess.item);//DA.3
            pManager.AddNumberParameter("Timestep", "Timestep", "Time division between the hours", GH_ParamAccess.item,1); //DA.4

        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddVectorParameter("Sun Vectors", "Vectors", "Sun Vectors for shadow or sunlight hour calculations", GH_ParamAccess.list);
        }

        /*
        protected override void AppendAdditionalComponentMenuItems(System.Windows.Forms.ToolStripDropDown menu)
        {
            base.AppendAdditionalComponentMenuItems(menu);
            Menu_AppendItem(menu, "Flip it",Menu_DoClick);
        }

        public static Mesh ms = new Mesh();

        private void Menu_DoClick(object sender, EventArgs e)
        {
            myBool = !myBool;
            RecordUndoEvent("Enabled Changed");
            ExpireSolution(true); //needed to get out from the class

        }
        
        public static bool myBool = false;
        */
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //0.0 Setting up the input parameters 

            Vector3d _North = Vector3d.YAxis;

            double lon = 0.0;
            double lat = 0.0;
            double tz = 0.0;
            double ts = 0.0;

            List<Vector3d> vecs = new List<Vector3d>();

            if (!DA.GetData(1, ref tz))
                return;

            if (!DA.GetData(2, ref lat))
                return;

            if (!DA.GetData(3, ref lon))
                return;

            try
            {
                DA.GetData(0, ref _North);
            }
            catch
            {
            }

            try
            {
                DA.GetData(4, ref ts);
            }
            catch
            {
            }

            //1.0 Main functions
            Vector3d North = new Vector3d(_North.X, _North.Y, 0.0);

            List<int> monthlist1 = new List<int>();
            monthlist1.Add(31);
            monthlist1.Add(28);
            monthlist1.Add(31);
            monthlist1.Add(30);
            monthlist1.Add(31);
            monthlist1.Add(30);
            monthlist1.Add(31);
            monthlist1.Add(31);
            monthlist1.Add(30);
            monthlist1.Add(31);
            monthlist1.Add(30);
            monthlist1.Add(31);

            List<int[]> setarray = new List<int[]>();

            for (int i = 0; i < 12; i++)
            {
                for (int j = 0; j < monthlist1[i]; j++)
                {
                    for (int h = 0; h < 24; h++)
                    {
                        for (int m = 0; m < ts; m++)
                        {
                            if (ts < 1) { ts = 1.0; }
                            int floor = Convert.ToInt32(Math.Floor(ts));
                            int incr = 60 / floor;

                            int[] toarr = new int[4] { h + 1, j + 1, i + 1, incr * m }; //HOUR (0), DAY(1), MONTH(2), MINS[3]
                            setarray.Add(toarr);
                        }
                    }
                }
            }

            var numbers = Enumerable.Range(0, setarray.Count());
            var _result = numbers.AsParallel().AsOrdered();
            foreach (var i in _result)
            {
                int[] place = setarray[i];


                double jday = getJD(place[2], place[1], 2017);

                double tl = getTimeLocal(place[0], place[3], 0);
                double total = jday + tl / 1440.0 - tz / 24.0;
                double T = calcTimeJulianCent(total);
                double _azimuth = calcAzEl(true, T, tl, lat, lon, tz);

                if(North.Length == 0) { North = Vector3d.YAxis; }

                Point3d p0 = new Point3d(0, 0, 0);
                Point3d p1 = new Point3d(0, 0, 0);
                p1.Transform(Transform.Translation(North));
                p1.Transform(Transform.Rotation(degToRad(-_azimuth), Vector3d.ZAxis, p0));

                Line l = new Line(p0, p1);

                Plane azPlane = new Plane(p0, l.Direction, Vector3d.ZAxis);
                p1.Transform(Transform.Rotation(degToRad(_altitude), azPlane.ZAxis, p0));

                Line d = new Line(p1, p0);
                Vector3d dir = d.Direction;
                dir.Unitize();
                //dir.Reverse();

                if (dir.Z <= 0)
                {
                    vecs.Add(dir);
                }
            }

            DA.SetDataList(0, vecs);

        }



        // Extra functions
        private static double alpha;
        private static double azRad;
        private static double _altitude;

        private static double refractionCorrection;

        private static double calcTimeJulianCent(double jd)
        {
            double T = (jd - 2451545.0) / 36525.0;
            return T;
        }

        private static double calcJDFromJulianCent(double t)
        {
            double JD = t * 36525.0 + 2451545.0;
            return JD;
        }

        private static bool isLeapYear(double yr)
        {
            return ((yr % 4 == 0 && yr % 100 != 0) || yr % 400 == 0);
        }

        private static List<int> monthlist = new List<int>();

        private static double getJD(int docmonth, int docday, double docyear)
        {
            monthlist.Add(31);
            monthlist.Add(28);
            monthlist.Add(31);
            monthlist.Add(30);
            monthlist.Add(31);
            monthlist.Add(30);
            monthlist.Add(31);
            monthlist.Add(31);
            monthlist.Add(30);
            monthlist.Add(31);
            monthlist.Add(30);
            monthlist.Add(31);

            if ((isLeapYear(docyear)) && (docmonth == 2))
            {
                if (docday > 29)
                {
                    docday = 29;
                    docday = docday - 1;
                }
            }
            else
            {
                if (docday > monthlist[docmonth - 1])
                {
                    docday = monthlist[docmonth - 1];
                }
            }

            if (docmonth <= 2)
            {
                docyear -= 1;
                docmonth += 12;
            }

            double A = Math.Floor(docyear / 100);
            double B = 2 - A + Math.Floor(A / 4);
            double JD = Math.Floor(365.25 * (docyear + 4716)) + Math.Floor(30.6001 * (docmonth + 1)) + docday + B - 1524.5;
            return JD;
        }

        private static double calcDoyFromJD(double jd)
        {
            double A;
            double z = Math.Floor(jd + 0.5);
            double f = (jd + 0.5) - z;
            if (z < 2299161)
            {
                A = z;
            }
            else
            {
                alpha = Math.Floor((z - 1867216.25) / 36524.25);
                A = z + 1 + alpha - Math.Floor(alpha / 4);
            }
            double B = A + 1524;
            double C = Math.Floor((B - 122.1) / 365.25);
            double D = Math.Floor(365.25 * C);
            double E = Math.Floor((B - D) / 30.6001);
            double day = B - D - Math.Floor(30.6001 * E) + f;
            double month = (E < 14) ? E - 1 : E - 13;
            double year = (month > 2) ? C - 4716 : C - 4715;

            int k = (isLeapYear(year) ? 1 : 2);
            double doy = Math.Floor((275 * month) / 9) - k * Math.Floor((month + 9) / 12) + day - 30;
            return doy;
        }

        private static double radToDeg(double angleRad)
        {
            return (180.0 * angleRad / Math.PI);
        }

        private static double degToRad(double angeDeg)
        {
            return (Math.PI * angeDeg / 180.0);
        }

        private static double calcGeomMeanLongSun(double t)
        {
            double L0 = 280.46646 + t * (36000.76983 + t * (0.0003032));
            while (L0 > 360.0)
            {
                L0 -= 360.0;
            }
            while (L0 < 0.0)
            {
                L0 += 360.0;
            }
            return L0;      // in degrees
        }

        private static double calcGeomMeanAnomalySun(double t)
        {
            double M = 357.52911 + t * (35999.05029 - 0.0001537 * t);
            return M;       // in degrees
        }

        private static double calcEccentricityEarthOrbit(double t)
        {
            double e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t);
            return e;       // unitless
        }

        private static double calcSunEqOfCenter(double t)
        {
            double m = calcGeomMeanAnomalySun(t);
            double mrad = degToRad(m);
            double sinm = Math.Sin(mrad);
            double sin2m = Math.Sin(mrad + mrad);
            double sin3m = Math.Sin(mrad + mrad + mrad);
            double C = sinm * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin2m * (0.019993 - 0.000101 * t) + sin3m * 0.000289;
            return C;       // in degrees
        }

        private static double calcSunTrueLong(double t)
        {
            double l0 = calcGeomMeanLongSun(t);
            double c = calcSunEqOfCenter(t);
            double O = l0 + c;
            return O;       // in degrees
        }

        private static double calcSunTrueAnomaly(double t)
        {
            double m = calcGeomMeanAnomalySun(t);
            double c = calcSunEqOfCenter(t);
            double v = m + c;
            return v;       // in degrees
        }

        private static double calcSunRadVector(double t)
        {
            double v = calcSunTrueAnomaly(t);
            double e = calcEccentricityEarthOrbit(t);
            double R = (1.000001018 * (1 - e * e)) / (1 + e * Math.Cos(degToRad(v)));
            return R;       // in AUs
        }

        private static double calcSunApparentLong(double t)
        {
            double o = calcSunTrueLong(t);
            double omega = 125.04 - 1934.136 * t;
            double lambda = o - 0.00569 - 0.00478 * Math.Sin(degToRad(omega));
            return lambda;      // in degrees
        }

        private static double calcMeanObliquityOfEcliptic(double t)
        {
            double seconds = 21.448 - t * (46.8150 + t * (0.00059 - t * (0.001813)));
            double e0 = 23.0 + (26.0 + (seconds / 60.0)) / 60.0;
            return e0;      // in degrees
        }

        private static double calcObliquityCorrection(double t)
        {
            double e0 = calcMeanObliquityOfEcliptic(t);
            double omega = 125.04 - 1934.136 * t;
            double e = e0 + 0.00256 * Math.Cos(degToRad(omega));
            return e;       // in degrees
        }

        private static double calcSunRtAscension(double t)
        {
            double e = calcObliquityCorrection(t);
            double lambda = calcSunApparentLong(t);
            double tananum = (Math.Cos(degToRad(e)) * Math.Sin(degToRad(lambda)));
            double tanadenom = (Math.Cos(degToRad(lambda)));
            double alp = radToDeg(Math.Atan2(tananum, tanadenom));
            return alp;       // in degrees
        }

        private static double calcSunDeclination(double t)
        {
            double e = calcObliquityCorrection(t);
            double lambda = calcSunApparentLong(t);

            double sint = Math.Sin(degToRad(e)) * Math.Sin(degToRad(lambda));
            double theta = radToDeg(Math.Asin(sint));
            return theta;       // in degrees
        }

        private static double calcEquationOfTime(double t)
        {
            double epsilon = calcObliquityCorrection(t);
            double l0 = calcGeomMeanLongSun(t);
            double e = calcEccentricityEarthOrbit(t);
            double m = calcGeomMeanAnomalySun(t);

            double y = Math.Tan(degToRad(epsilon) / 2.0);
            y *= y;

            double sin2l0 = Math.Sin(2.0 * degToRad(l0));
            double sinm = Math.Sin(degToRad(m));
            double cos2l0 = Math.Cos(2.0 * degToRad(l0));
            double sin4l0 = Math.Sin(4.0 * degToRad(l0));
            double sin2m = Math.Sin(2.0 * degToRad(m));

            double Etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0 - 0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m;
            return radToDeg(Etime) * 4.0;   // in minutes of time
        }

        private static double calcHourAngleSunrise(double lat, double solarDec)
        {
            double latRad = degToRad(lat);
            double sdRad = degToRad(solarDec);
            double HAarg = (Math.Cos(degToRad(90.833)) / (Math.Cos(latRad) * Math.Cos(sdRad)) - Math.Tan(latRad) * Math.Tan(sdRad));
            double HA = Math.Acos(HAarg);
            return HA;      // in radians (for sunset, use -HA)
        }

        private static double getTimeLocal(int dochr, int docmn, int docsc)
        {

            double mins = dochr * 60 + docmn + docsc / 60.0;
            return mins;
        }

        private static double calcAzEl(bool output, double T, double localtime, double latitude, double longitude, double zone)
        {
            double eqTime = calcEquationOfTime(T);
            double theta = calcSunDeclination(T);
            double azimuth;

            double solarTimeFix = eqTime + 4.0 * longitude - 60.0 * zone;
            double earthRadVec = calcSunRadVector(T);
            double trueSolarTime = localtime + solarTimeFix;

            while (trueSolarTime > 1440)
            {
                trueSolarTime -= 1440;
            }


            double hourAngle = trueSolarTime / 4.0 - 180.0;
            if (hourAngle < -180)
            {
                hourAngle += 360.0;
            }
            double haRad = degToRad(hourAngle);

            double csz = Math.Sin(degToRad(latitude)) * Math.Sin(degToRad(theta)) + Math.Cos(degToRad(latitude)) * Math.Cos(degToRad(theta)) * Math.Cos(haRad);
            if (csz > 1.0)
            {
                csz = 1.0;
            }

            else if (csz < -1.0)
            {
                csz = -1.0;
            }

            double zenith = radToDeg(Math.Acos(csz));
            double azDenom = (Math.Cos(degToRad(latitude)) * Math.Sin(degToRad(zenith)));

            if (Math.Abs(azDenom) > 0.001)
            {
                azRad = ((Math.Sin(degToRad(latitude)) * Math.Cos(degToRad(zenith))) - Math.Sin(degToRad(theta))) / azDenom;
                if (Math.Abs(azRad) > 1.0)
                {
                    if (azRad < 0)
                    {
                        azRad = -1.0;
                    }
                    else
                    {
                        azRad = 1.0;
                    }
                }
                azimuth = 180.0 - radToDeg(Math.Acos(azRad));
                if (hourAngle > 0.0)
                {
                    azimuth = -azimuth;
                }
            }
            else
            {
                if (latitude > 0.0)
                {
                    azimuth = 180.0;
                }
                else
                {
                    azimuth = 0.0;
                }
            }
            if (azimuth < 0.0)
            {
                azimuth += 360.0;
            }
            double exoatmElevation = 90.0 - zenith;

            // Atmospheric Refraction correction

            if (exoatmElevation > 85.0)
            {
                refractionCorrection = 0.0;
            }
            else
            {
                var te = Math.Tan(degToRad(exoatmElevation));
                if (exoatmElevation > 5.0)
                {
                    refractionCorrection = 58.1 / te - 0.07 / (te * te * te) + 0.000086 / (te * te * te * te * te);
                }
                else if (exoatmElevation > -0.575)
                {
                    refractionCorrection = 1735.0 + exoatmElevation * (-518.2 + exoatmElevation * (103.4 + exoatmElevation * (-12.79 + exoatmElevation * 0.711)));
                }
                else
                {
                    refractionCorrection = -20.774 / te;
                }
                refractionCorrection = refractionCorrection / 3600.0;
            }

            var solarZen = zenith - refractionCorrection;

            _altitude = Math.Floor((90.0 - solarZen) * 100 + 0.5) / 100.0;
            azimuth = Math.Floor(azimuth * 100 + 0.5) / 100.0;
            return azimuth;
        }



        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                return Properties.Resources.solarcalc;
                //return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("3F9BA182-03F6-491E-A3E1-EF359B722F9A"); }
        }

    }
    
    public class MarsSolarPath : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public MarsSolarPath()
          : base("Martian Solar Path", "Martian Solar",
              "Give solar vectors in Mars, given any coordinates",
              "BIG", "Analysis")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Latitude", "Lat", "North coordinate of the location on Mars", GH_ParamAccess.item);
            pManager.AddNumberParameter("Longitude", "Lon", "East coordinate positions of the location on Mars", GH_ParamAccess.item);
            pManager.AddVectorParameter("North", "North", "North location of the Rhino", GH_ParamAccess.item,Vector3d.YAxis);

            pManager.AddNumberParameter("Day", "Day", "Optional day input (int 1 to 669) to show sun movement for that day", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("Hour", "Hour", "Optional Martian hour input (int 1 to 24) to show sun movement for that hour", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("Radius", "Radius", "Optional Radius for the sunpath, default is 100 units", GH_ParamAccess.item, 100.0);

            pManager.AddPointParameter("Center", "Center", "Optional center location for the sunpath, default is origin", GH_ParamAccess.item, Point3d.Origin);
            pManager.AddBooleanParameter("NightTime", "Night", "Set true to include the sun position during the night time", GH_ParamAccess.item, false);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            /*
             Param01 = Vectors for sunlight hours analysis
             Param02 = Curves solely for nice and beautz diagrams
             Param03 = Sunposition for rhino light creation
             Param04 = Coordinated Local Mars Time 
             Param05 = Martian Sol Azimuth
             Param06 = Martian Sol Altitude 
             */

            pManager.AddVectorParameter("SunVectors", "Vectors", "Sun incoming angle to the mars surface", GH_ParamAccess.list); 
            pManager.AddCurveParameter("SunCurves", "Curves", "Sun curves and compass", GH_ParamAccess.list);
            pManager.AddPointParameter("SunPositions", "Points", "Sun positions of given period", GH_ParamAccess.list);
            pManager.AddTextParameter("LMST", "LMST", "Coordinated local mean solar time", GH_ParamAccess.list);
            pManager.AddNumberParameter("Azimuth", "Azi", "Azimuth of martian solar (in degrees)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Altitude", "Alt", "Altitude of martian solar (in degrees)", GH_ParamAccess.list);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            //--------------------- Output containers 

            List<Vector3d> vecs = new List<Vector3d>();
            List<Curve> crv_diag = new List<Curve>();
            List<Point3d> pts_sun = new List<Point3d>();

            List<string> LMST = new List<string>();
            List<double> _azi = new List<double>();
            List<double> _alt = new List<double>();

            //--------------------- Initial settings

            bool night = false;
            double lat = 0.0;
            double lon = 0.0;
            double dayz = 0;
            double hourz = 0;
            double radius = 100.0;
            Point3d center = Point3d.Origin;
            Vector3d _North = Vector3d.YAxis;

            List<string> _LMST = new List<string>();
            List<double> _azi0 = new List<double>();
            List<double> _alt0 = new List<double>();

            if (!DA.GetData(0, ref lat))
            {
                return;
            }

            if (!DA.GetData(1, ref lon))
            {
                return;
            }

            if (!DA.GetData(2, ref _North))
            {
                return;
            }

            if (!DA.GetData(3, ref dayz))
            {
                return;
            }

            if (!DA.GetData(4, ref hourz))
            {
                return;
            }

            if (!DA.GetData(5, ref radius))
            {
                return;
            }

            if(!DA.GetData(6,ref center))
            {
                return;
            }

            dayz = dayz % 669;
            hourz = hourz % 24;

            if (!DA.GetData(7, ref night))
            {
                return;
            }

            Vector3d North = new Vector3d(_North.X, _North.Y, 0);
            North.Unitize();
            if (North.Length == 0) { North = Vector3d.YAxis; }
            

            //--------------------- Produce Azimuth and Altitude

            if (dayz == 0)
            {
                for (int i = 1; i < 670; i++)
                {
                    if (hourz == 0)
                    {
                        for (int j = 1; j < 25; j++)
                        {
                            double jdttoff = j2000_ott_from_Mars_Solar_Date(i) + (double)j / 24;
                            _alt.Add(solar_elevation(lon, lat, jdttoff));
                            _azi.Add(solar_azimuth(lon, lat, jdttoff));
                            _LMST.Add(GenTimeSpanFromHours( Local_Mean_Solar_Time(lon,jdttoff)));

                        }
                    }
                    else
                    {
                        double jdttoff = j2000_ott_from_Mars_Solar_Date(i) + hourz / 24;
                        _alt.Add(solar_elevation(lon, lat, jdttoff));
                        _azi.Add(solar_azimuth(lon, lat, jdttoff));
                        _LMST.Add(GenTimeSpanFromHours(Local_Mean_Solar_Time(lon, jdttoff)));

                    }
                }
            }
            else
            {
                if (hourz == 0)
                {
                    for (int j = 1; j < 25; j++)
                    {
                        double jdttoff = j2000_ott_from_Mars_Solar_Date(dayz) + (double)j / 24;
                        _alt.Add(solar_elevation(lon, lat, jdttoff));
                        _azi.Add(solar_azimuth(lon, lat, jdttoff));
                        _LMST.Add(GenTimeSpanFromHours(Local_Mean_Solar_Time(lon, jdttoff)));

                    }
                }
                else
                {
                    double jdttoff = j2000_ott_from_Mars_Solar_Date(dayz) + hourz / 24;
                    _alt.Add(solar_elevation(lon, lat, jdttoff));
                    _azi.Add(solar_azimuth(lon, lat, jdttoff));
                    _LMST.Add(GenTimeSpanFromHours(Local_Mean_Solar_Time(lon, jdttoff)));

                }
            }

            //--------------------- Produce Vectors

            for (int i = 0; i<_alt.Count();i++)
            {
                double _azimuth = _azi[i];
                double _altitude = _alt[i];

                Vector3d dir = altaztovecs(_altitude, _azimuth, North);

                if (!night)
                {
                    if (dir.Z <= 0)
                    {
                        vecs.Add(dir);
                        LMST.Add(_LMST[i]);
                        _azi0.Add(_azimuth);
                        _alt0.Add(_altitude);
                    }
                }
                else
                {
                    vecs.Add(dir);
                }


            }

            if (night)
            {
                _azi0.Clear();
                _alt0.Clear();
                LMST.Clear();

                LMST = _LMST;
                _azi0 = _azi;
                _alt0 = _alt;
            }


            //--------------------- Produce SunPosition
            foreach(Vector3d v in vecs)
            {
                double _rad = -radius;
                Line sunline = new Line(center, v, _rad);
                pts_sun.Add(sunline.PointAt(1.0));
            }

            //--------------------- Produce SunCurve

            //Main Compass
            Circle comp_crl = new Circle(center, 1.1*radius);
            crv_diag.Add(comp_crl.ToNurbsCurve());

            //Secondary Compass
            comp_crl.Transform(Transform.Scale(center, 1/1.1));
            crv_diag.Add(comp_crl.ToNurbsCurve());

            //Lines 
            Point3d pNorth0 = new Point3d(center.X, center.Y, center.Z);
            pNorth0.Transform(Transform.Translation(1 * radius * North));

            Point3d pNorth1 = new Point3d(pNorth0.X, pNorth0.Y, pNorth0.Z);
            pNorth1.Transform(Transform.Translation(0.22 * radius * North));

            Point3d pNorth2 = new Point3d(pNorth0.X, pNorth0.Y, pNorth0.Z);
            pNorth2.Transform(Transform.Translation(0.1 * radius * North));

            Line lnnorth = new Line(pNorth0, pNorth1);
            crv_diag.Add(lnnorth.ToNurbsCurve()); //north

            Line lnnorthsmall = new Line(pNorth0, pNorth2);
            for (int i = 1; i < 16; i++) //all other 8 directions
            {
                lnnorthsmall.Transform(Transform.Rotation(degToRad(i*22.5), Vector3d.ZAxis, center));
                crv_diag.Add(lnnorthsmall.ToNurbsCurve()); 
            }

            //8 days curve 

            int[] daystodiag = new int[8] { 1, 83, 167, 251, 335, 418, 502, 585 };
            foreach(int it in daystodiag)
            {
                List<Point3d> pts_inside = new List<Point3d>();
                for(int i = 1; i < 4; i++)
                {
                    double jdttoff = j2000_ott_from_Mars_Solar_Date(it) + (double) i / 3;
                    Vector3d dir = altaztovecs(solar_elevation(lon, lat, jdttoff), solar_azimuth(lon, lat, jdttoff), North);
                    Line sunline = new Line(center, dir, -radius);
                    pts_inside.Add(sunline.PointAt(1.0));
                }

                Circle cex = new Circle(pts_inside[0], pts_inside[1], pts_inside[2]);

                if (!night) //if noon only
                {
                    NurbsCurve split = cex.ToNurbsCurve();
                    var ccx = Rhino.Geometry.Intersect.Intersection.CurvePlane(split, Plane.WorldXY, 0.02);
                    if (ccx.Count() < 2)
                    {
                        crv_diag.Add(cex.ToNurbsCurve());
                    }
                    else 
                    {
                        double[] dbsplit = new double[2] { ccx[0].ParameterA, ccx[1].ParameterA };
                        Curve[] splits = new Curve[0];
                        splits = split.Split(dbsplit);

                        foreach (Curve c in splits)
                        {
                            Curve cs0 = c;
                            cs0.Domain = new Interval(0, 1);
                            if (cs0.PointAt(0.5).Z > 0) { crv_diag.Add(cs0); }
                        }

                    }
                    Message = "Only daytime";

                }
                else //if night
                {
                    Message = "Night time also";
                    crv_diag.Add(cex.ToNurbsCurve());
                }
            }

            /*
             Param01 = Vectors for sunlight hours analysis
             Param02 = Curves solely for nice and beautz diagrams
             Param03 = Sunposition for rhino light creation
             Param04 = Coordinated Local Mars Time 
             Param05 = Martian Sol Azimuth
             Param06 = Martian Sol Altitude 
             */

            DA.SetDataList(0, vecs);//done
            DA.SetDataList(1, crv_diag); //must check the nighttime 
            DA.SetDataList(2, pts_sun);//done

            DA.SetDataList(3, LMST); //done
            DA.SetDataList(4, _azi0); //done
            DA.SetDataList(5, _alt0); //done

        }

        // X. Extra additional useful functions
        private static Vector3d altaztovecs (double alt, double az, Vector3d North)
        {

            Point3d p0 = new Point3d(0, 0, 0);
            Point3d p1 = new Point3d(0, 0, 0);
            Point3d p2 = new Point3d(0, 0, 1);

            p1.Transform(Transform.Translation(North));
            p1.Transform(Transform.Rotation(degToRad(-az), Vector3d.ZAxis, p0));

            Line l = new Line(p0, p1);

            // Plane azPlane = new Plane(p0, l.Direction, Vector3d.ZAxis);
            Plane azPlane = new Plane(p0, p1, p2);
            p1.Transform(Transform.Rotation(degToRad(alt), azPlane.ZAxis, p0));

            Line d = new Line(p1, p0);
            Vector3d dir = d.Direction;
            dir.Unitize();

            return dir;
        }
        private static double west_to_east(double west)
        {
            //Convert from west longitude to east longitude,
            //or vice versa. 
            double east = 360 - west;
            //This function used to convert between two coordinate systems at the same time,
            //which seems to be wrong, so I've removed that option.
            return east % 360;
        }

        private static string GenTimeSpanFromHours(double hours)
        {
            // Create a TimeSpan object and TimeSpan string from 
            // a number of hours.
            TimeSpan interval = TimeSpan.FromHours(hours);
            string timeInterval = interval.ToString(@"hh\:mm\:ss");

            return timeInterval;

        }

        private static double east_to_west(double east)
        {
            return west_to_east(east);
        }

        private static double j2000_epoch = 2451545.0;

        private static double radToDeg(double angleRad)
        {
            return (180.0 * angleRad / Math.PI);
        }

        private static double degToRad(double angeDeg)
        {
            return (Math.PI * angeDeg / 180.0);
        }

        private static double mills()
        {
            DateTime today = DateTime.Today;
            DateTime epoch = new DateTime(1970, 1, 1, 0, 0, 0);
            TimeSpan ts = today - epoch;
            return ts.Milliseconds;
        }

        private static double julian(double mill)
        {
            return 2440587.5 + (mill / 8.64e7);
        }

        private static double utc_to_tt_offset(double jday)
        {
            return utc_to_tt_offset_math(jday);
        }

        //input julian day to return the universal time constant
        private static double utc_to_tt_offset_math(double jday)
        {
            double jday_min = 2441317.5;
            double offset_min = 32.184;

            double jday_np = jday;

            double[] jday_vals = new double[27] { -2441317.5, 0,    182,    366,
                       731,   1096,   1461,   1827,
                       2192,   2557,   2922,   3469,
                       3834,   4199,   4930,   5844,
                       6575,   6940,   7487,   7852,
                       8217,   8766,   9313,   9862,
                       12419,  13515, 14792 };

            double[] offset_vals = new double[27] {-32.184,10.0, 11.0, 12.0, 13.0,
                    14.0, 15.0, 16.0, 17.0, 18.0,
                    19.0, 20.0, 21.0, 22.0, 23.0,
                    24.0, 25.0, 26.0, 27.0, 28.0,
                    29.0, 30.0, 31.0, 32.0, 33.0,
                    34.0, 35.0 };



            if (jday_np <= jday_min + jday_vals[0])
            {
                return offset_min + offset_vals[0];
            }
            else if (jday_np >= jday_min + jday_vals[jday_vals.Count()-1])
            {
                return offset_min + offset_vals[offset_vals.Count() - 1];
            }
            else
            {
                int i = 0;
                for (i = 0; i < offset_vals.Count(); i++)
                {

                    if (jday_min + jday_vals[i] <= jday_np && jday_min + jday_vals[i + 1] > jday_np)
                    {
                        break;
                    }
                }
                return offset_min + offset_vals[i];
            }
        }

        private static double julian_tt(double jday_utc)
        {

            double jdtt = jday_utc + utc_to_tt_offset(jday_utc) / 86400.0;
            return jdtt;

        }

        private static double j2000_offset_tt(double jday_tt)
        {
            return jday_tt - 2451545.0;
        }

        private static double Mars_Mean_Anomaly(double j2000_ott)
        {
            double M = 19.3870 + 0.52402075 * j2000_ott;
            return M % 360;
        }

        private static double FMS_Angle(double j2000_ott)
        {
            double alpha_fms = 270.3863 + 0.52403840 * j2000_ott;
            return alpha_fms % 360;
        }

        private static double alpha_perturbs(double j2000_ott)
        {
            double[] array_A = new double[7] { 0.0071, 0.0057, 0.0039, 0.0037, 0.0021, 0.0020, 0.0018 };
            double[] array_tau = new double[7] { 2.2353, 2.7543, 1.1177, 15.7866, 2.1354, 2.4694, 32.8493 };
            double[] array_phi = new double[7] { 49.409, 168.173, 191.837, 21.736, 15.704, 95.528, 49.095 };

            double pbs = 0;


            for(int i = 0; i< array_A.Count(); i++)
            {
                pbs = pbs + array_A[i] * Math.Cos(((0.985626 * j2000_ott / array_tau[i]) + array_phi[i]) * Math.PI / 180.0);
            }

            return pbs;
        }

        private static double equation_of_center(double j2000_ott)
        {
            double M = Mars_Mean_Anomaly(j2000_ott) * Math.PI / 180;
            double pbs = alpha_perturbs(j2000_ott);

            double val = (10.691 + 3.0e-7 * j2000_ott) * Math.Sin(M) +
                0.6230 * Math.Sin(2 * M) +
                0.0500 * Math.Sin(3 * M) +
                0.0050 * Math.Sin(4 * M) +
                0.0005 * Math.Sin(5 * M) + pbs;
            return val;
        }

        private static double Mars_Ls(double j2000_ott)
        {
            double alpha = FMS_Angle(j2000_ott);
            double v_m = equation_of_center(j2000_ott);
            double ls = (alpha + v_m);
            ls = ls % 360;
            return ls;
        }

        private static double equation_of_time(double j2000_ott)
        {
            double ls = Mars_Ls(j2000_ott) * Math.PI / 180;
            double EOT = 2.861 * Math.Sin(2 * ls) -
                0.071 * Math.Sin(4 * ls) +
                0.002 * Math.Sin(6 * ls) - equation_of_center(j2000_ott);

            return EOT;

        }

        private static double j2000_from_Mars_Solar_Date(double msd)
        {
            double j2000_ott = ((msd + 0.00096 - 44796.0) * 1.027491252) + 4.5;
            return j2000_ott;
        }


        private static double j2000_ott_from_Mars_Solar_Date(double msd)
        {
            double j2000 = j2000_from_Mars_Solar_Date(msd);
            double j2000_ott = julian_tt(j2000 + j2000_epoch);
            return j2000_ott - j2000_epoch;

        }

        private static double Mars_Solar_Date(double j2000_ott)
        {
            double MSD = (((j2000_ott - 4.5) / 1.027491252) + 44796.0 - 0.00096);
            return MSD;
        }

        private static double Clancy_Year(double j2000_ott)
        {
            double ref1955_4_11_11am = -16336.0416;
            double year = Math.Floor(1 + (j2000_ott - ref1955_4_11_11am) / (686.978));
            return year;
        }

        private static double[] Mars_Year(double j2000_ott, bool return_length)
        {
            double[] jday_vals = new double[79] { -16336.044076, -15649.093471, -14962.0892946, -14275.0960023, -13588.1458658, -12901.1772635, -12214.2082215, -11527.2637345, -10840.2842249, -10153.2828749, -9466.3114025, -8779.3356111, -8092.3607738, -7405.4236452, -6718.4615347, -6031.4574604, -5344.4876509, -4657.5318339, -3970.5474528, -3283.5848372, -2596.6329362, -1909.6426682, -1222.6617049, -535.7040268, 151.2736522, 838.2369682, 1525.1834712, 2212.1799182, 2899.1848518, 3586.1403058, 4273.1024234, 4960.0765368, 5647.0207838, 6333.986502, 7020.9875066, 7707.9629132, 8394.9318782, 9081.9102062, 9768.8526533, 10455.8028354, 11142.8050514, 11829.7873254, 12516.7417734, 13203.725449, 13890.6991502, 14577.6484912, 15264.6324865, 15951.6217969, 16638.5798914, 17325.5517216, 18012.5209097, 18699.4628887, 19386.4443201, 20073.4534421, 20760.4152811, 21447.3696661, 22134.3466251, 22821.2966642, 23508.2529432, 24195.2539572, 24882.2400506, 25569.2081296, 26256.1902459, 26943.1429481, 27630.0847446, 28317.0793316, 29004.0710936, 29691.0238241, 30377.9991486, 31064.9784277, 31751.9249377, 32438.896907, 33125.8902412, 33812.8520242, 34499.8183442, 35186.7944595, 35873.740573, 36560.7112423, 37247.7247318 };
            double[] year_vals = new double[79] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79 };
            double[] year_length = new double[79] { 686.95252, 686.950605, 687.0041764, 686.9932923, 686.9501365, 686.9686023, 686.969042, 686.944487, 686.9795096, 687.00135, 686.9714724, 686.9757914, 686.9748373, 686.9371286, 686.9621105, 687.0040743, 686.9698095, 686.955817, 686.9843811, 686.9626156, 686.951901, 686.990268, 686.9809633, 686.9576781, 686.977679, 686.963316, 686.946503, 686.996447, 687.0049336, 686.955454, 686.9621176, 686.9741134, 686.944247, 686.9657182, 687.0010046, 686.9754066, 686.968965, 686.978328, 686.9424471, 686.9501821, 687.002216, 686.982274, 686.954448, 686.9836756, 686.9737012, 686.949341, 686.9839953, 686.9893104, 686.9580945, 686.9718302, 686.9691881, 686.941979, 686.9814314, 687.009122, 686.961839, 686.954385, 686.976959, 686.9500391, 686.956279, 687.001014, 686.9860934, 686.968079, 686.9821163, 686.9527022, 686.9417965, 686.994587, 686.991762, 686.9527305, 686.9753245, 686.9792791, 686.94651, 686.9719693, 686.9933342, 686.961783, 686.96632, 686.9761153, 686.9461135, 686.9706693, 687.0134895 };

            return Mars_Year_math(j2000_ott, jday_vals, year_vals, year_length, return_length);
        }

        
        private static double[] Mars_Year_math(double j2k_math, double[] jday_vals, double[] year_vals, double[] year_length, bool return_length)
        {
            int j = 0;
            int right = jday_vals.Count() - 1;
            if (j2k_math < jday_vals[0])
            {
                double[] dbs = new double[1] { Math.Floor(1 + (j2k_math - jday_vals[0]) / year_length[0]) };
                return dbs;
            }
            else if (j2k_math >= jday_vals[right])
            {
                double[] dbs = new double[1] { Math.Floor(1 + (j2k_math - jday_vals[right]) / year_length[right]) };
                return dbs;
            }
            else
            {
                for (int i = 0; i < year_vals.Count() - 1; i++)
                {
                    if ((jday_vals[i] <= j2k_math) && (jday_vals[i + 1] > j2k_math))
                    {
                        j = i;
                        break;
                    }
                }
                
            }

            double y = year_vals[j];
            double l = year_length[j];

            if (return_length)
            {
                double[] dbs = new double[2] { y, l };
                return dbs;
            }
            else
            {
                double[] dbs = new double[1] { y };
                return dbs;
            }
        }
        
        private static double Coordinated_Mars_Time(double j2000_ott)
        {
            double MTC = 24 * (((j2000_ott - 4.5) / 1.027491252) + 44796.0 - 0.00096);
            MTC = MTC % 24;
            return MTC;
        }

        private static double Local_Mean_Solar_Time(double longitude, double j2000_ott)
        {
            double MTC = Coordinated_Mars_Time(j2000_ott);
            double LMST = MTC - longitude * (24 / 360);
            LMST = LMST % 24;
            return LMST;
        }

        private static double Local_True_Solar_Time(double longitude, double j2000_ott)
        {
            double eot = equation_of_time(j2000_ott);
            double lmst = Local_Mean_Solar_Time(longitude, j2000_ott);
            double ltst = lmst + eot * (24 / 360);
            ltst = ltst % 24;
            return ltst;
        }

        private static double subsolar_longitude(double j2000_ott)
        {
            double MTC = Coordinated_Mars_Time(j2000_ott);
            double EOT = equation_of_time(j2000_ott) * 24 / 360;
            double subsol = (MTC + EOT) * (360 / 24) + 180;
            return subsol % 360;
        }

        private static double solar_declination(double ls)
        {
            double ls1 = ls * Math.PI / 180;
            double dec = Math.Asin(0.42565 * Math.Sin(ls1)) + 0.25 * (Math.PI / 180) * Math.Sin(ls1);
            dec = dec * 180 / Math.PI;
            return dec;
        }

        private static double heliocentric_distance(double j2000_ott)
        {
            double M = Mars_Mean_Anomaly(j2000_ott) * Math.PI / 180;
            double rm = 1.523679 *(1.00436 - 0.09309 * Math.Cos(M)
             - 0.004336 * Math.Cos(2 * M)
             - 0.00031 * Math.Cos(3 * M)
             - 0.00003 * Math.Cos(4 * M));
            return rm;
        }

        private static double heliocentric_longitude(double j2000_ott)
        {
            double ls = Mars_Ls(j2000_ott);
            double im = ls + 85.061 - 0.015 * Math.Sin((71 + 2 * ls) * Math.PI / 180) - 
        5.5e-6 * j2000_ott;
            return im % 360;
        }

        private static double heliocentric_latitude(double j2000_ott)
        {
            double ls = Mars_Ls(j2000_ott);
            double bm = -(1.8497 - 2.23e-5 * j2000_ott) 
        * Math.Sin((ls - 144.50 + 2.57e-6 * j2000_ott) * Math.PI / 180);
            return bm;
        }

        private static double hourangle(double longitude, double j2000_ott)
        {
            double subsol = subsolar_longitude(j2000_ott) * Math.PI / 180;
            double hourangle = longitude * Math.PI / 180 - subsol;
            return hourangle;
        }

        private static double solar_zenith(double longitude, double latitude, double j2000_ott)
        {
            double ha = hourangle(longitude, j2000_ott);
            double ls = Mars_Ls(j2000_ott);
            double dec = solar_declination(ls) * Math.PI / 180;

            double cosZ = Math.Sin(dec) * Math.Sin(latitude * Math.PI/ 180) + Math.Cos(dec) * Math.Cos(latitude * Math.PI / 180) * Math.Cos(ha);

            return Math.Acos(cosZ) * 180/ Math.PI;
        }

        private static double solar_elevation(double longitude, double latitude, double j2000_ott)
        {
            return 90 - solar_zenith(longitude, latitude, j2000_ott);
        }

        private static double solar_azimuth(double longitude, double latitude, double j2000_ott)
        {
            double ha = hourangle(longitude, j2000_ott);
            double ls = Mars_Ls(j2000_ott);
            double dec = solar_declination(ls) * Math.PI / 180;
            double denom = (Math.Cos(degToRad( latitude)) * Math.Tan(dec)
                 - Math.Sin(degToRad(latitude)) * Math.Cos(ha));

            double num = Math.Sin(ha);

            return (360 + Math.Atan2(num, denom) * 180/ Math.PI) % 360;
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                return Properties.Resources.mars;
                //return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("90193628-8E62-4EBB-9475-FA106E671105"); }
        }
    }
    
    public class VisualizeSun : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public VisualizeSun()
          : base("VisualizeSun", "VizSun",
              "Show the sun points from given sun vectors",
              "BIG", "Analysis")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddVectorParameter("Sun Vectors", "SunVectors", "Vectors from the sun components", GH_ParamAccess.item);
            pManager.AddPointParameter("Center", "Center", "Center point for the sunpath", GH_ParamAccess.item, Point3d.Origin);
            pManager.AddNumberParameter("Length", "Length", "Distance from the center point", GH_ParamAccess.item, 1000);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Sun", "Sun", "Sun location", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d origin = Point3d.Origin;
            double length = 10000;
            Vector3d vec = new Vector3d(0, 0, 1);

            if (!DA.GetData(0, ref vec))
            {
                return;
            }

            if (!DA.GetData(1, ref origin))
            {
                return;
            }

            if (!DA.GetData(2, ref length))
            {
                return;
            }

            Line l = new Line(origin, vec, -length);
            Point3d sunis = l.PointAt(1.0);

            DA.SetData(0, sunis);

        }


        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                return Properties.Resources._2_2_sun_high_quality_png;
                //return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("DFE4ACE9-DA84-49E7-A6CF-D7EB5DF053B2"); }
        }
    }
    
}
