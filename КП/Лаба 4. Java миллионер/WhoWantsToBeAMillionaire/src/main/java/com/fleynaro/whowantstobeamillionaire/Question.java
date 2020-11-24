/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.fleynaro.whowantstobeamillionaire;

/**
 *
 * @author user
 */
public class Question
{
    public Question(String[] s)
    {
        Text = s[0];
        for (int i=0; i<4; i++)
            Answers[i]= s[i+1];
        RightAnswer=s[5];
        Level = Integer.parseInt(s[6]);
    }
    
    public String Text;
    String[] Answers = new String[4];
    public String RightAnswer;
    public int Level;
}
