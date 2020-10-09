using System;
using System.Windows.Forms;
using System.Drawing;

namespace ControlsLab1
{
    public partial class ColorSelector : UserControl
    {
        public Color Color
        {
            get { return color; }
            set { ChangeColor(value, true); }
        }

        public event EventHandler ColorChanged;

        private Color color;
        public ColorSelector()
        {
            InitializeComponent();
            color = Color.FromArgb(255, 0, 0, 0);
            colorSelectorInput_R.ColorChanged += ColorSelectorInput_SomeColorComponentChanged;
            colorSelectorInput_B.ColorChanged += ColorSelectorInput_SomeColorComponentChanged;
            colorSelectorInput_G.ColorChanged += ColorSelectorInput_SomeColorComponentChanged;
        }

        private void ColorSelectorInput_SomeColorComponentChanged(byte value)
        {
            var color = Color.FromArgb(255, colorSelectorInput_R.ColorValue, colorSelectorInput_G.ColorValue, colorSelectorInput_B.ColorValue);
            ChangeColor(color);
        }

        private void radioButton_Dec_CheckedChanged(object sender, EventArgs e)
        {
            colorSelectorInput_R.Mode =
            colorSelectorInput_G.Mode =
            colorSelectorInput_B.Mode = ColorSelectorInput.EnumMode.Dec;
        }

        private void radioButton_Hex_CheckedChanged(object sender, EventArgs e)
        {
            colorSelectorInput_R.Mode =
            colorSelectorInput_G.Mode =
            colorSelectorInput_B.Mode = ColorSelectorInput.EnumMode.Hex;
        }

        protected override void OnPaint(PaintEventArgs e)
        {
            base.OnPaint(e);
            SolidBrush brush = new SolidBrush(color);
            e.Graphics.FillRectangle(brush, GetColorViewRectangle());
        }

        private void ChangeColor(Color color, bool updateInputs = false)
        {
            this.color = color;
            if (updateInputs)
            {
                colorSelectorInput_R.ColorValue = color.R;
                colorSelectorInput_G.ColorValue = color.G;
                colorSelectorInput_B.ColorValue = color.B;
            }
            ColorChanged?.Invoke(this, new EventArgs());
            Invalidate();
        }

        private Rectangle GetColorViewRectangle()
        {
            int X = colorSelectorInput_R.Location.X + colorSelectorInput_R.Size.Width + 10;
            int Y = colorSelectorInput_R.Location.Y;
            return new Rectangle(new Point(X, Y), new Size(100, 100));
        }
    }
}
