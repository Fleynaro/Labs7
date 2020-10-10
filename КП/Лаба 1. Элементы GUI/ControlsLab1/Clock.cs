using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra;

namespace ControlsLab1
{
    public partial class Clock : Control
    {
        public class Appearance
        {
            public Color BackgroundColor = Color.Gray;
            public Color ClockHandColor = Color.Black;
            public Color ClockFaceColor = Color.DarkGray;
        }

        public Appearance App
        {
            get { return app; }

            set { app = value; Invalidate(); }
        }

        Appearance app = new Appearance();
        Vector<float> clockDefSize;
        Vector<float> clockDefPos;
        Matrix<float> clockHourHandMatrix;
        Matrix<float> clockMinHandMatrix;
        Matrix<float> clockSecHandMatrix;
        Matrix<float> clockHourLabelMatrix;
        Matrix<float> clockMinLabelMatrix;

        public Clock()
        {
            InitializeComponent();
            clockDefSize = Vector<float>.Build.DenseOfArray(new float[] { 200.0f, 200.0f });
            clockDefPos = Vector<float>.Build.DenseOfArray(new float[] { 0.0f, 0.0f });
            clockHourHandMatrix = CreateRectangleMatrix(new SizeF(5.0f, 50.0f));
            clockMinHandMatrix = CreateRectangleMatrix(new SizeF(3.0f, 70.0f));
            clockSecHandMatrix = CreateRectangleMatrix(new SizeF(2.0f, 70.0f));
            clockHourLabelMatrix = CreateRectangleMatrix(new SizeF(3.0f, 20.0f));
            clockMinLabelMatrix = CreateRectangleMatrix(new SizeF(2.0f, 10.0f));
        }

        protected override void OnPaint(PaintEventArgs pe)
        {
            base.OnPaint(pe);
            var scaleMatrix = GetScaleMatrix();
            
            //draw circle
            var brush = new SolidBrush(app.BackgroundColor);
            pe.Graphics.FillEllipse(brush, new RectangleF(
                new PointF(clockDefPos[0], clockDefPos[1]),
                new SizeF(Width, Height))
            );

            
            var clockCenterPos = clockDefPos + clockDefSize / 2;
            var centerTranslMatrix = Matrix<float>.Build.DenseOfColumnVectors(
                clockCenterPos, clockCenterPos, clockCenterPos, clockCenterPos);

            Matrix<float> matrix;
            var PI2 = (float)Math.PI * 2;

            //draw labels
            brush = new SolidBrush(app.ClockFaceColor);
            var labelOffset = Vector<float>.Build.DenseOfArray(new float[] { 0, clockDefSize[1] * 0.48f });
            var labelOffsetMatrix = Matrix<float>.Build.DenseOfColumnVectors(
                labelOffset, labelOffset, labelOffset, labelOffset);
            for (int hourIdx = 0; hourIdx < 12; hourIdx++)
            {
                matrix = (CreateRotationMatrix(hourIdx * PI2 / 12) * (clockHourLabelMatrix + labelOffsetMatrix)) + centerTranslMatrix;
                DrawRectangle(pe, brush, scaleMatrix * matrix);
            }
            for (int minIdx = 0; minIdx < 60; minIdx++)
            {
                if (minIdx % 10 == 0)
                    continue;
                matrix = (CreateRotationMatrix(minIdx * PI2 / 60) * (clockMinLabelMatrix + labelOffsetMatrix)) + centerTranslMatrix;
                DrawRectangle(pe, brush, scaleMatrix * matrix);
            }


            //get components of time
            var dateTime = GetDateTime();
            var hour = dateTime.Hour % 12;
            var minute = dateTime.Minute;
            var second = dateTime.Second;

            //draw hour hands
            brush = new SolidBrush(app.ClockHandColor);
            //hour hand
            matrix = (CreateRotationMatrix(hour * PI2 / 12) * clockHourHandMatrix) + centerTranslMatrix;
            DrawRectangle(pe, brush, scaleMatrix * matrix);
            //minute hand
            matrix = (CreateRotationMatrix(minute * PI2 / 60) * clockMinHandMatrix) + centerTranslMatrix;
            DrawRectangle(pe, brush, scaleMatrix * matrix);
            //sec hand
            matrix = (CreateRotationMatrix(second * PI2 / 60) * clockSecHandMatrix) + centerTranslMatrix;
            DrawRectangle(pe, brush, scaleMatrix * matrix);
        }

        protected DateTime GetDateTime()
        {
            return DateTime.Now;
        }

        private void DrawRectangle(PaintEventArgs pe, Brush brush, Matrix<float> matrix)
        {
            var vertexes = new PointF[] {
                new PointF(matrix[0, 0], matrix[1, 0]),
                new PointF(matrix[0, 1], matrix[1, 1]),
                new PointF(matrix[0, 2], matrix[1, 2]),
                new PointF(matrix[0, 3], matrix[1, 3])
            };
            pe.Graphics.FillPolygon(brush, vertexes);
        }

        private Matrix<float> CreateRectangleMatrix(SizeF size)
        {
            var widthVec = Vector<float>.Build.DenseOfArray(new float[] { size.Width, 0 });
            var heightVec = Vector<float>.Build.DenseOfArray(new float[] { 0, -size.Height });
            var leftDownPos = -widthVec / 2;
            var leftUpPos = leftDownPos + heightVec;
            var rightDownPos = -leftDownPos;
            var rightUpPos = rightDownPos + heightVec;
            return Matrix<float>.Build.DenseOfColumnVectors(
                leftDownPos,
                leftUpPos,
                rightUpPos,
                rightDownPos
            );
        }

        private Matrix<float> CreateRotationMatrix(float alpha)
        {
            return Matrix<float>.Build.DenseOfArray(new float[,] {
                { (float)Math.Cos(alpha), -(float)Math.Sin(alpha) },
                { (float)Math.Sin(alpha), (float)Math.Cos(alpha) }
            });
        }

        private Matrix<float> GetScaleMatrix()
        {
            return Matrix<float>.Build.DenseOfArray(new float[,] {
                { Width / clockDefSize[0], 0 },
                { 0, Height / clockDefSize[1] }
            });
        }
    }
}
