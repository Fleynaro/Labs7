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
            this.clock1 = new ControlsLab1.Clock();
            this.colorSelector1 = new ControlsLab1.ColorSelector();
            this.SuspendLayout();
            // 
            // clock1
            // 
            this.clock1.Location = new System.Drawing.Point(315, 37);
            this.clock1.Name = "clock1";
            this.clock1.Size = new System.Drawing.Size(250, 231);
            this.clock1.TabIndex = 1;
            this.clock1.Text = "clock1";
            // 
            // colorSelector1
            // 
            this.colorSelector1.Color = System.Drawing.Color.FromArgb(((int)(((byte)(0)))), ((int)(((byte)(0)))), ((int)(((byte)(0)))));
            this.colorSelector1.Location = new System.Drawing.Point(26, 80);
            this.colorSelector1.Name = "colorSelector1";
            this.colorSelector1.Size = new System.Drawing.Size(255, 113);
            this.colorSelector1.TabIndex = 0;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(594, 305);
            this.Controls.Add(this.clock1);
            this.Controls.Add(this.colorSelector1);
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);

        }

        #endregion

        private ControlsLab1.ColorSelector colorSelector1;
        private ControlsLab1.Clock clock1;
    }
}

