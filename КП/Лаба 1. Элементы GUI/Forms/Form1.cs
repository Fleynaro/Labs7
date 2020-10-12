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
            
        }

        private void ColorSelector1_ColorChanged(object sender, EventArgs e)
        {
            clockWPF1.BorderColor = colorSelector1.Color;
        }
    }
}
