namespace Forms
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.colorSelector1 = new ControlsLab1.ColorSelector();
            this.elementHost1 = new System.Windows.Forms.Integration.ElementHost();
            this.clockWPF1 = new ControlsLab1.ClockWPF();
            this.SuspendLayout();
            // 
            // colorSelector1
            // 
            this.colorSelector1.Color = System.Drawing.Color.FromArgb(((int)(((byte)(0)))), ((int)(((byte)(0)))), ((int)(((byte)(0)))));
            this.colorSelector1.Location = new System.Drawing.Point(26, 80);
            this.colorSelector1.Name = "colorSelector1";
            this.colorSelector1.Size = new System.Drawing.Size(255, 113);
            this.colorSelector1.TabIndex = 0;
            // 
            // elementHost1
            // 
            this.elementHost1.Location = new System.Drawing.Point(314, 28);
            this.elementHost1.Name = "elementHost1";
            this.elementHost1.Size = new System.Drawing.Size(460, 297);
            this.elementHost1.TabIndex = 1;
            this.elementHost1.Text = "elementHost1";
            this.elementHost1.Child = this.clockWPF1;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(837, 410);
            this.Controls.Add(this.elementHost1);
            this.Controls.Add(this.colorSelector1);
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);

        }

        #endregion

        private ControlsLab1.ColorSelector colorSelector1;
        private System.Windows.Forms.Integration.ElementHost elementHost1;
        private ControlsLab1.ClockWPF clockWPF1;
    }
}

