using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace ComForms
{
    public partial class Form1 : Form
    {
        ATM atm;
        public Form1()
        {
            InitializeComponent();
            atm = new ATM();
            ICreateAtm crAtm = (ICreateAtm)atm;
            crAtm.SetLocationAddr("г. Пермь, ул. Культуры, д. 10");
            crAtm.SetMoney(1000);
        }

        private void button2_Click(object sender, EventArgs e)
        {
            IStats stat = (IStats)atm;
            stat.DisplayStats();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            IBalance balance = (IBalance)atm;
            if( radioButton1.Checked)
            {
                balance.PayMoney(0, (int)numericUpDown1.Value);
                MessageBox.Show("Деньги успешно положены на карту!");
            }

            if (radioButton2.Checked)
            {
                balance.WithdrawMoney(0, (int)numericUpDown1.Value);
                MessageBox.Show("Деньги успешно сняты с карты!");
            }

            numericUpDown1.Value = 0;
        }
    }
}
