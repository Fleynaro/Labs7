namespace ControlsLab1
{
    partial class ColorSelector
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

        #region Component Designer generated code

        /// <summary> 
        /// Required method for Designer support - do not modify 
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.label1 = new System.Windows.Forms.Label();
            this.radioButton_Dec = new System.Windows.Forms.RadioButton();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.radioButton_Hex = new System.Windows.Forms.RadioButton();
            this.colorSelectorInput_R = new ControlsLab1.ColorSelectorInput(this.components);
            this.colorSelectorInput_G = new ControlsLab1.ColorSelectorInput(this.components);
            this.colorSelectorInput_B = new ControlsLab1.ColorSelectorInput(this.components);
            this.SuspendLayout();
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(16, 21);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(52, 13);
            this.label1.TabIndex = 1;
            this.label1.Text = "Красный";
            // 
            // radioButton_Dec
            // 
            this.radioButton_Dec.AutoSize = true;
            this.radioButton_Dec.Location = new System.Drawing.Point(19, 97);
            this.radioButton_Dec.Name = "radioButton_Dec";
            this.radioButton_Dec.Size = new System.Drawing.Size(45, 17);
            this.radioButton_Dec.TabIndex = 2;
            this.radioButton_Dec.TabStop = true;
            this.radioButton_Dec.Text = "Dec";
            this.radioButton_Dec.UseVisualStyleBackColor = true;
            this.radioButton_Dec.CheckedChanged += new System.EventHandler(this.radioButton_Dec_CheckedChanged);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(16, 48);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(52, 13);
            this.label2.TabIndex = 4;
            this.label2.Text = "Зелёный";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(16, 74);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(38, 13);
            this.label3.TabIndex = 6;
            this.label3.Text = "Синий";
            // 
            // radioButton_Hex
            // 
            this.radioButton_Hex.AutoSize = true;
            this.radioButton_Hex.Location = new System.Drawing.Point(76, 97);
            this.radioButton_Hex.Name = "radioButton_Hex";
            this.radioButton_Hex.Size = new System.Drawing.Size(44, 17);
            this.radioButton_Hex.TabIndex = 7;
            this.radioButton_Hex.TabStop = true;
            this.radioButton_Hex.Text = "Hex";
            this.radioButton_Hex.UseVisualStyleBackColor = true;
            this.radioButton_Hex.CheckedChanged += new System.EventHandler(this.radioButton_Hex_CheckedChanged);
            // 
            // colorSelectorInput_R
            // 
            this.colorSelectorInput_R.Location = new System.Drawing.Point(76, 19);
            this.colorSelectorInput_R.Name = "colorSelectorInput_R";
            this.colorSelectorInput_R.Size = new System.Drawing.Size(45, 20);
            this.colorSelectorInput_R.TabIndex = 8;
            // 
            // colorSelectorInput_G
            // 
            this.colorSelectorInput_G.Location = new System.Drawing.Point(75, 45);
            this.colorSelectorInput_G.Name = "colorSelectorInput_G";
            this.colorSelectorInput_G.Size = new System.Drawing.Size(45, 20);
            this.colorSelectorInput_G.TabIndex = 9;
            // 
            // colorSelectorInput_B
            // 
            this.colorSelectorInput_B.Location = new System.Drawing.Point(75, 71);
            this.colorSelectorInput_B.Name = "colorSelectorInput_B";
            this.colorSelectorInput_B.Size = new System.Drawing.Size(45, 20);
            this.colorSelectorInput_B.TabIndex = 10;
            // 
            // ColorSelector
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.Controls.Add(this.colorSelectorInput_B);
            this.Controls.Add(this.colorSelectorInput_G);
            this.Controls.Add(this.colorSelectorInput_R);
            this.Controls.Add(this.radioButton_Hex);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.radioButton_Dec);
            this.Controls.Add(this.label1);
            this.Name = "ColorSelector";
            this.Size = new System.Drawing.Size(250, 120);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.RadioButton radioButton_Dec;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.RadioButton radioButton_Hex;
        private ColorSelectorInput colorSelectorInput_R;
        private ColorSelectorInput colorSelectorInput_G;
        private ColorSelectorInput colorSelectorInput_B;
    }
}
