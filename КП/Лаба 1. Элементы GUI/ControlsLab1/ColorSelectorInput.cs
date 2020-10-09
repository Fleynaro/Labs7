using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Windows.Forms;
using System.Globalization;

namespace ControlsLab1
{
    public partial class ColorSelectorInput : TextBox
    {
        public byte ColorValue
        {
            get { return colorValue; }

            set { ChangeColorValue(value); }
        }

        public EnumMode Mode
        {
            set { ChangeMode(value); }
        }

        public enum EnumMode {
            Dec,
            Hex
        }

        public delegate void ColorChangeEventHandler(byte value);
        public event ColorChangeEventHandler ColorChanged;

        private byte colorValue;
        private bool isValid;
        private EnumMode mode;

        public ColorSelectorInput()
        {
            InitializeComponent();
            Text = "0";
            colorValue = 0;
            isValid = true;
            mode = EnumMode.Dec;
        }

        public ColorSelectorInput(IContainer container)
        {
            container.Add(this);

            InitializeComponent();
        }

        protected override void OnTextChanged(EventArgs e)
        {
            base.OnTextChanged(e);

            int value;
            var format = (mode == EnumMode.Dec) ? NumberStyles.Integer : NumberStyles.HexNumber;
            isValid = int.TryParse(Text, format, CultureInfo.CurrentCulture, out value);
            if (isValid)
            {
                if (value < 0) value = 0;
                if (value > 255) value = 255;
                colorValue = (byte)value;
                ColorChanged?.Invoke(colorValue);
            }

            if (!isValid)
            {
                ForeColor = Color.Red;
            }
            else
            {
                ForeColor = Color.Black;
            }
        }

        protected override void OnKeyPress(KeyPressEventArgs e)
        {
            if (mode == EnumMode.Dec && !char.IsDigit(e.KeyChar) && !char.IsControl((e.KeyChar)))
                e.Handled = true;
            base.OnKeyPress(e);
        }
        
        private void ChangeMode(EnumMode mode)
        {
            this.mode = mode;
            UpdateText();
        }

        private void ChangeColorValue(byte value)
        {
            this.colorValue = value;
            UpdateText();
        }

        private void UpdateText()
        {
            if (mode == EnumMode.Hex)
            {
                Text = colorValue.ToString("X");
            }
            else
            {
                Text = colorValue.ToString();
            }
        }
    }
}
