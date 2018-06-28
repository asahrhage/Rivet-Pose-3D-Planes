using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenCvSharp;
using MathNet.Spatial;
using MathNet.Spatial.Euclidean;
using MathNet.Spatial.Units;

namespace Rivet3DPlanes
{
    class Program
    {

        static Point3D _cameraCenter1;
        static Point3D _cameraCenter2;


        static Mat _R;
        static Mat _T;
        static Mat _RInverse;

        static Mat _cameraMatrixLeft;

        static double _fxLeft;
        static double _fyLeft;
        static double _fLeft;
        static double _cxLeft;
        static double _cyLeft;




        static Mat _cameraMatrixRight;

        static double _fxRight;
        static double _fyRight;
        static double _fRight;
        static double _cxRight;
        static double _cyRight;



        static List<Point3d> _rivetImagePointsPlaneLeft;
        static List<Point3d> _rivetImagePointsPlaneRight;

        static List<Point3D> _rivetImagePointsPlaneLeftCam1COS;
        static List<Point3D> _rivetImagePointsPlaneRightCam1COS;

        static List<Point> _rivetBoundingBoxPointsImageLeft;
        static List<Point> _rivetBoundingBoxPointsImageRight;

        static List<Point> _rivetAxisPointsImageLeft;
        static List<Point> _rivetAxisPointsImageRight;

        static double _rivetRadius = 2;

        static Point3D _rivetCenter1;
        static Point3D _rivetCenter2;



        static void Main(string[] args)
        {

            ReadCameraMatrix();
            ReadRivets();
            
            CalculateRivetPose3D(_cameraCenter1, _rivetImagePointsPlaneLeftCam1COS[0], _rivetImagePointsPlaneLeftCam1COS[1], _cameraCenter2, _rivetImagePointsPlaneRightCam1COS[0], _rivetImagePointsPlaneRightCam1COS[1]);

            Console.WriteLine(_rivetCenter1.X.ToString() + ";" + _rivetCenter1.Y.ToString() + ";" + _rivetCenter1.Z.ToString());
            Console.WriteLine(_rivetCenter2.X.ToString() + ";" + _rivetCenter2.Y.ToString() + ";" + _rivetCenter2.Z.ToString());

            Console.ReadKey();



        }






        static void ReadCameraMatrix()
        {




            string filepathToXML = @"C:\Users\Arne\Documents\Masterarbeit\IMAGES\stereoSetup3\calibration\maps\maps.xml";
            using (FileStorage fs = new FileStorage(filepathToXML, FileStorage.Mode.Read))
            {
                _cameraMatrixLeft = fs["cam1Mat"].ReadMat();
                _cameraMatrixRight = fs["cam2Mat"].ReadMat();

            }

            //double[,] cameraMatrixArrayLeft = new double[3, 3];
            ////cameraMatrixArray[0, 0] = 3085.4f;
            ////cameraMatrixArray[0, 1] = 0f;
            ////cameraMatrixArray[0, 2] = 573.32f;
            ////cameraMatrixArray[1, 0] = 0;
            ////cameraMatrixArray[1, 1] = 3076.5f;
            ////cameraMatrixArray[1, 2] = 632.41f;
            ////cameraMatrixArray[2, 0] = 0;
            ////cameraMatrixArray[2, 1] = 0;
            ////cameraMatrixArray[2, 2] = 1;

            //cameraMatrixArrayLeft[0, 0] = 3085.38539652417;
            //cameraMatrixArrayLeft[0, 1] = 0;
            //cameraMatrixArrayLeft[0, 2] = 573.216605993667;
            //cameraMatrixArrayLeft[1, 0] = 0;
            //cameraMatrixArrayLeft[1, 1] = 3076.50773957751;
            //cameraMatrixArrayLeft[1, 2] = 632.406561211840;
            //cameraMatrixArrayLeft[2, 0] = 0;
            //cameraMatrixArrayLeft[2, 1] = 0;
            //cameraMatrixArrayLeft[2, 2] = 1;

            //_cameraMatrixLeft = new Mat(3, 3, MatType.CV_64FC1, cameraMatrixArrayLeft);


            //_fxLeft = cameraMatrixArrayLeft[0, 0];
            //_fyLeft = cameraMatrixArrayLeft[1, 1];
            //_fLeft = (_fxLeft + _fyLeft) / 2;


            //_cxLeft = cameraMatrixArrayLeft[0, 2];
            //_cyLeft = cameraMatrixArrayLeft[1, 2];

            _fxLeft = _cameraMatrixLeft.GetArray(0, 0)[0];
            _fyLeft = _cameraMatrixLeft.GetArray(1, 1)[0];
            _fLeft = (_fxLeft + _fyLeft) / 2;


            _cxLeft = _cameraMatrixLeft.GetArray(1, 2)[0];
            _cyLeft = _cameraMatrixLeft.GetArray(0,2)[0];

            //double[,] cameraMatrixArrayRight = new double[3, 3];
            ////cameraMatrixArray[0, 0] = 3085.4f;
            ////cameraMatrixArray[0, 1] = 0f;
            ////cameraMatrixArray[0, 2] = 573.32f;
            ////cameraMatrixArray[1, 0] = 0;
            ////cameraMatrixArray[1, 1] = 3076.5f;
            ////cameraMatrixArray[1, 2] = 632.41f;
            ////cameraMatrixArray[2, 0] = 0;
            ////cameraMatrixArray[2, 1] = 0;
            ////cameraMatrixArray[2, 2] = 1;

            //cameraMatrixArrayRight[0, 0] = 3068.69478897440;
            //cameraMatrixArrayRight[0, 1] = 0;
            //cameraMatrixArrayRight[0, 2] = 496.134171106576;
            //cameraMatrixArrayRight[1, 0] = 0;
            //cameraMatrixArrayRight[1, 1] = 3064.134473839508;
            //cameraMatrixArrayRight[1, 2] = 633.010075465329;
            //cameraMatrixArrayRight[2, 0] = 0;
            //cameraMatrixArrayRight[2, 1] = 0;
            //cameraMatrixArrayRight[2, 2] = 1;

            //_cameraMatrixRight = new Mat(3, 3, MatType.CV_64FC1, cameraMatrixArrayRight);


            //_fxRight = cameraMatrixArrayRight[0, 0];
            //_fyRight = cameraMatrixArrayRight[1, 1];
            //_fRight = (_fxRight + _fyRight) / 2;


            //_cxRight = cameraMatrixArrayRight[0, 2];
            //_cyRight = cameraMatrixArrayRight[1, 2];

            _fxRight = _cameraMatrixRight.GetArray(0, 0)[0];
            _fyRight = _cameraMatrixRight.GetArray(1, 1)[0];
            _fRight = (_fxRight + _fyRight) / 2;


            _cxRight = _cameraMatrixRight.GetArray(1, 2)[0];
            //_cyRight = 548;
            _cyRight = _cameraMatrixRight.GetArray(0, 2)[0];



            //string filepathToXML = @"C:\Users\Arne\Documents\Masterarbeit\IMAGES\stereoSetup3\calibration\maps\maps.xml";
            using (FileStorage fs = new FileStorage(filepathToXML, FileStorage.Mode.Read))
            {
                _R = fs["R"].ReadMat();
                _T = fs["T"].ReadMat();

            }


            double[,] RArray = new double[3, 3];
            //RArray[0, 0] = 0.933630108366741;
            //RArray[0, 1] = 0.0488847095357635;
            //RArray[0, 2] = -0.354887455293523;
            //RArray[1, 0] = -0.0638646141839124;
            //RArray[1, 1] = 0.997488942207857;
            //RArray[1, 2] = -0.0306124358421566;
            //RArray[2, 0] = 0.352499832349248;
            //RArray[2, 1] = 0.0512454422037138;
            //RArray[2, 2] = 0.934407712322141;

            //RArray[0, 0] = 0.933630108366741;
            //RArray[0, 1] = 0.0488847095357635;
            //RArray[0, 2] = -0.354887455293523;
            //RArray[1, 0] = -0.0638646141839124;
            //RArray[1, 1] = 0.997488942207857;
            //RArray[1, 2] = -0.0306124358421566;
            //RArray[2, 0] = 0.352499832349248;
            //RArray[2, 1] = 0.0512454422037138;
            //RArray[2, 2] = 0.934407712322141;

            //_R = new Mat(3, 3, MatType.CV_64FC1, RArray);

            //double[] TArray = new double[3];
            //TArray[0] = -41.7547848045858;
            //TArray[1] = -2.49332048778849;
            //TArray[2] = 7.00137268354101;

            //_T = new Mat(3, 1, MatType.CV_64FC3, TArray);


            _cameraCenter1 = new Point3D(0, 0, 0);
            _cameraCenter2 = new Point3D(-1*_T.GetArray(0,0)[0], -1 * _T.GetArray(1, 0)[0], -1 * _T.GetArray(2, 0)[0]);


            _RInverse = _R.Inv();


        }



        static void ReadRivets()
        {
            _rivetBoundingBoxPointsImageLeft = new List<Point>();
            _rivetBoundingBoxPointsImageRight = new List<Point>();
            _rivetImagePointsPlaneLeft = new List<Point3d>();
            _rivetImagePointsPlaneRight = new List<Point3d>();



            //Add points for middle axis 

            // __x__
            //|     |
            //|     |
            //|__x__| 

            //_rivetBoundingBoxPointsImageLeft.Add(new Point(319, 598));
            //_rivetBoundingBoxPointsImageLeft.Add(new Point(312, 467));
            //_rivetBoundingBoxPointsImageLeft.Add(new Point(488, 457));
            //_rivetBoundingBoxPointsImageLeft.Add(new Point(495, 588));

            //_rivetBoundingBoxPointsImageRight.Add(new Point(175,674));
            //_rivetBoundingBoxPointsImageRight.Add(new Point(175, 546));
            //_rivetBoundingBoxPointsImageRight.Add(new Point(380, 546));
            //_rivetBoundingBoxPointsImageRight.Add(new Point(380, 674));


            _rivetBoundingBoxPointsImageLeft.Add(new Point(319, 598));
            _rivetBoundingBoxPointsImageLeft.Add(new Point(312, 467));
            _rivetBoundingBoxPointsImageLeft.Add(new Point(488, 457));
            _rivetBoundingBoxPointsImageLeft.Add(new Point(495, 588));

            _rivetBoundingBoxPointsImageRight.Add(new Point(175, 674));
            _rivetBoundingBoxPointsImageRight.Add(new Point(175, 546));
            _rivetBoundingBoxPointsImageRight.Add(new Point(380, 546));
            _rivetBoundingBoxPointsImageRight.Add(new Point(380, 674));


            _rivetAxisPointsImageLeft = GetMiddleAxisFromBoundingBox(_rivetBoundingBoxPointsImageLeft);
            _rivetAxisPointsImageRight = GetMiddleAxisFromBoundingBox(_rivetBoundingBoxPointsImageRight);

            //_rivetImageRight.Add(new Point(0, 0));
            //_rivetImageRight.Add(new Point(0, 0));


            ReprojectPointsImageImagePlane(_rivetAxisPointsImageLeft, _rivetImagePointsPlaneLeft, _fxLeft, _fyLeft, _fLeft, _cxLeft, _cyLeft);
            ReprojectPointsImageImagePlane(_rivetAxisPointsImageRight, _rivetImagePointsPlaneRight, _fxRight, _fyRight, _fRight, _cxRight, _cyRight);



            _rivetImagePointsPlaneLeftCam1COS = new List<Point3D>();
            _rivetImagePointsPlaneLeftCam1COS.Add(new Point3D(_rivetImagePointsPlaneLeft[0].X, _rivetImagePointsPlaneLeft[0].Y, _rivetImagePointsPlaneLeft[0].Z));
            _rivetImagePointsPlaneLeftCam1COS.Add(new Point3D(_rivetImagePointsPlaneLeft[1].X, _rivetImagePointsPlaneLeft[1].Y, _rivetImagePointsPlaneLeft[1].Z));

            _rivetImagePointsPlaneRightCam1COS =  TransformPointsToCam1COS(_rivetImagePointsPlaneRight, _R, _T);




        }

        static void ReprojectPointsImageImagePlane(List<Point> imageCoordinates, List<Point3d> imagePlaneCoordinates, double fx, double fy, double f, double cx, double cy)
        {

            for (int i = 0; i < imageCoordinates.Count(); i++)
            {

                
                double xTempLeft = (imageCoordinates[i].X - cx) * 5.3 / 1000;
                double yTempLeft = (imageCoordinates[i].Y - cy) * 5.3 / 1000;
                double zTempLeft = f*5.3/1000;
                imagePlaneCoordinates.Add(new Point3d(xTempLeft, yTempLeft, zTempLeft));

                




            }



            
        }


        static void CalculateRivetPose3D(Point3D cameraCenter1, Point3D pointCam1_1, Point3D pointCam1_2,Point3D cameraCenter2,  Point3D pointCam2_1, Point3D pointCam2_2)
        {




            //Plane planeCam1 = new Plane(new Point3D(cameraCenter1.X, cameraCenter1.Y, cameraCenter1.Z), new Point3D(pointCam1_1.X, pointCam1_1.Y, pointCam1_1.Z), new Point3D(pointCam1_2.X, pointCam1_2.Y, pointCam1_2.Z));
            Plane planeCam2 = new Plane(new Point3D(cameraCenter2.X, cameraCenter2.Y, cameraCenter2.Z), new Point3D(pointCam2_1.X, pointCam2_1.Y, pointCam2_1.Z), new Point3D(pointCam2_2.X, pointCam2_2.Y, pointCam2_2.Z));

            Vector3D VectorCam1CenterToPoint1_1 = new Vector3D(pointCam1_1.X - cameraCenter1.X, pointCam1_1.Y - cameraCenter1.Y, pointCam1_1.Z - cameraCenter1.Z);
            Vector3D VectorCam1CenterToPoint1_2 = new Vector3D(pointCam1_2.X - cameraCenter1.X, pointCam1_2.Y - cameraCenter1.Y, pointCam1_2.Z - cameraCenter1.Z);

            //Ray3D rivetAxis3D =  planeCam1.IntersectionWith(planeCam2);

            Ray3D RayCam1CenterToPoint1_1 = new Ray3D(cameraCenter1, VectorCam1CenterToPoint1_1.Normalize());
            Ray3D RayCam1CenterToPoint1_2 = new Ray3D(cameraCenter1, VectorCam1CenterToPoint1_2.Normalize());


            Point3D pointRivetCylinder1 = planeCam2.IntersectionWith(RayCam1CenterToPoint1_1);
            Point3D pointRivetCylinder2 = planeCam2.IntersectionWith(RayCam1CenterToPoint1_2);

            //var a = rivetAxis3D.LineTo(pointRivetCylinder1);
            //var b = rivetAxis3D.LineTo(pointRivetCylinder2);


            //Angle alpha1 = RayCam1CenterToPoint1_1.Direction.AngleTo(rivetAxis3D.Direction);
            //Angle alpha2 = RayCam1CenterToPoint1_2.Direction.AngleTo(rivetAxis3D.Direction);

            UnitVector3D intersection1To2 = new Vector3D(pointRivetCylinder2.X - pointRivetCylinder1.X, pointRivetCylinder2.Y - pointRivetCylinder1.Y, pointRivetCylinder2.Z - pointRivetCylinder1.Z).Normalize();

            Angle alpha1 = RayCam1CenterToPoint1_1.Direction.AngleTo(intersection1To2);
            Angle alpha2 = RayCam1CenterToPoint1_2.Direction.AngleTo(intersection1To2);

            if (pointRivetCylinder1.DistanceTo(new Point3D(0,0,0)) < pointRivetCylinder2.DistanceTo(new Point3D(0, 0, 0)))
            {

                double distance1 = Math.Abs(_rivetRadius / Math.Tan(alpha1.Radians));
                Vector3D distanceAlongRivetAxis1 = distance1 * intersection1To2;
                _rivetCenter1 = pointRivetCylinder1 + distanceAlongRivetAxis1;
                double distance2 = Math.Abs(_rivetRadius / Math.Tan(alpha2.Radians));
                Vector3D distanceAlongRivetAxis2 = distance2 * intersection1To2;
                _rivetCenter2 = pointRivetCylinder2 - distanceAlongRivetAxis2;




            }

            else
            {
                double distance1 = _rivetRadius / Math.Tan(alpha1.Radians);
                Vector3D distanceAlongRivetAxis1 = distance1 * intersection1To2;
                _rivetCenter1 = pointRivetCylinder1 - distanceAlongRivetAxis1;
                double distance2 = _rivetRadius / Math.Tan(alpha2.Radians);
                Vector3D distanceAlongRivetAxis2 = distance2 * intersection1To2;
                _rivetCenter2 = pointRivetCylinder2 + distanceAlongRivetAxis2;
            }



            //MathNet.Spatial.Units.Angle plane1AxisAngle = planeCam1.Normal.AngleTo(rivetAxis3D.Direction);
            double rivetLength = Math.Sqrt(Math.Pow(_rivetCenter1.X - _rivetCenter2.X, 2) + Math.Pow(_rivetCenter1.Y - _rivetCenter2.Y, 2) + Math.Pow(_rivetCenter1.Z - _rivetCenter2.Z, 2));

            //Console.WriteLine("distance between points: " + rivetLength.ToString());


            //var length1 = rivetAxis3D.LineTo(pointRivetCylinder1).Length;
            //var length2 = rivetAxis3D.LineTo(pointRivetCylinder2).Length;

            //var collinear = rivetAxis3D.Direction.IsParallelTo(intersection1To2);


        }



        static List<Point3D> TransformPointsToCam1COS(List<Point3d> pointsCam2COS, Mat R, Mat T)
        {
            List<Point3D> pointsCam1COS = new List<Point3D>();



            foreach (var pointCam2COS in pointsCam2COS)
            {
                double[] pointArray = new double[3] { pointCam2COS.X, pointCam2COS.Y, pointCam2COS.Z };
                Mat pointMat = new Mat(3, 1, MatType.CV_64FC1, pointArray);


                //Mat TNeg = -1 * T;
                Mat pointCam1COSMat = _RInverse*pointMat;
                //pointCam1COSMat = pointCam1COSMat + TNeg;
                Point3D tempPoint = new Point3D(pointCam1COSMat.GetArray(0, 0)[0]-T.GetArray(0,0)[0], pointCam1COSMat.GetArray(1, 0)[0] - T.GetArray(1, 0)[0], pointCam1COSMat.GetArray(2, 0)[0] - T.GetArray(2, 0)[0]);
                pointsCam1COS.Add(tempPoint);

            }




            double[] pointAxisArray = new double[3] { 0, 0, 1 };
            Mat pointAxisMat = new Mat(3, 1, MatType.CV_64FC1, pointAxisArray);

            Mat pointAxisCam1COSMAt = _RInverse * pointAxisMat;
            Point3D pointAxisCAM1COS = new Point3D(pointAxisCam1COSMAt.GetArray(0, 0)[0], pointAxisCam1COSMAt.GetArray(1,0)[0], pointAxisCam1COSMAt.GetArray(2,0)[0]);



            return pointsCam1COS;



        }


        static List<Point> GetMiddleAxisFromBoundingBox(List<Point> boundingBox)
        {

            List<Point> pointsMiddleAxis = new List<Point>();
            List<double> distances = new List<double>();
            for (int i = 1; i < boundingBox.Count(); i++)
            {
                double dist = Math.Sqrt(Math.Pow(boundingBox[0].X - boundingBox[i].X, 2) + Math.Pow(boundingBox[0].Y - boundingBox[i].Y, 2));
                distances.Add(dist);
            }

            for (int i = 0; i < boundingBox.Count(); i++)
            {

                for (int j = i+1; j < boundingBox.Count(); j++)
                {
                    double dist = Math.Sqrt(Math.Pow(boundingBox[i].X - boundingBox[j].X, 2) + Math.Pow(boundingBox[i].Y - boundingBox[j].Y, 2));


                    if (dist==distances.Min())
                    {

                        Point tempMiddle = new Point((int)((boundingBox[i].X + boundingBox[j].X) / 2), (int)((boundingBox[i].Y + boundingBox[j].Y) / 2));
                        pointsMiddleAxis.Add(tempMiddle);


                    }




                }

            }




            return pointsMiddleAxis;
        }



    }


}
