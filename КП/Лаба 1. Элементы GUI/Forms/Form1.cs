using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Forms
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            colorSelector1.ColorChanged += ColorSelector1_ColorChanged;
            colorSelector1.Color = clock1.App.BackgroundColor;

            var timer = new Timer();
            timer.Interval = 1000;
            timer.Tick += Timer_Tick;
            timer.Start();
        }

        private void Timer_Tick(object sender, EventArgs e)
        {
            clock1.Invalidate();
        }

        private void ColorSelector1_ColorChanged(object sender, EventArgs e)
        {
            var app = clock1.App;
            app.BackgroundColor = colorSelector1.Color;
            clock1.App = app;
        }
    }
}
