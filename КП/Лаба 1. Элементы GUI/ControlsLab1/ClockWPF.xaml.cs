using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Threading;

namespace ControlsLab1
{
    /// <summary>
    /// Interaction logic for ClockWPF.xaml
    /// </summary>
    public partial class ClockWPF : UserControl
    {
        public int TimeZone = 0;
        public System.Drawing.Color BorderColor
        {
            set { ClockBorder.Fill = new SolidColorBrush(Color.FromArgb(value.A, value.R, value.G, value.B)); }
        }

        private DispatcherTimer timer = new DispatcherTimer();

        public ClockWPF()
        {
            InitializeComponent();
            timer.Interval = new TimeSpan(0, 0, 1);
            timer.Tick += OnUpdate;
            timer.Start();
            Update();
        }

        protected DateTime GetDateTime()
        {
            return DateTime.Now;
        }

        private void OnUpdate(object sender, object e)
        {
            Update();
        }

        public void Update()
        {
            var time = GetDateTime();
            var secHandAngle = time.Second * 6.0;
            var minHandAngle = time.Minute * 6.0 + secHandAngle / 60.0;
            var hourHandAngle = 30.0 * ((time.Hour + TimeZone) % 12) + minHandAngle / 12.0;
            SecHand.RenderTransform = new RotateTransform(secHandAngle);
            MinHand.RenderTransform = new RotateTransform(minHandAngle);
            HourHand.RenderTransform = new RotateTransform(hourHandAngle);
        }
    }
}
