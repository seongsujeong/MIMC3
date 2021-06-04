#ifndef _MIMC2_MISC_
#include "MIMC_misc.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

char isLeapYear(int yint)
{
    char flag_leap=0;
    if(yint%4==0)
    {
        if(yint%100==0)
        {
            if(yint%400==0)
            {
                flag_leap=1;
            }
        }
        else
        {
            flag_leap=1;
        }
    }
    return flag_leap;
}

double get_datenum(char *string_date)
{
    //float out;
    double offset_date;
    int yint,mint,dint;
    double Hint,Mint,Sint;
    int vec_accumulate_date[12];
    int cnt_year;
    char flag_leap=0;
    char string_y[5],string_m[3],string_d[3],string_H[3],string_M[3],string_S[3];

    string_y[4]='\0';
    string_m[2]='\0';
    string_d[2]='\0';
    string_H[2]='\0';
    string_M[2]='\0';
    string_S[2]='\0';
    
    //parse the input string
    string_y[0]=string_date[0];
    string_y[1]=string_date[1];
    string_y[2]=string_date[2];
    string_y[3]=string_date[3];

    string_m[0]=string_date[4];
    string_m[1]=string_date[5];

    string_d[0]=string_date[6];
    string_d[1]=string_date[7];

    string_H[0]=string_date[8];
    string_H[1]=string_date[9];

    string_M[0]=string_date[10];
    string_M[1]=string_date[11];

    string_S[0]=string_date[12];
    string_S[1]=string_date[13];

    //convert the string to integer
    yint=atoi(string_y);
    mint=atoi(string_m);
    dint=atoi(string_d);
    Hint=atof(string_H);
    Mint=atof(string_M);
    Sint=atof(string_S);

    offset_date=0.0;
    for(cnt_year=0;cnt_year<yint;cnt_year++)
    {
        flag_leap=isLeapYear(cnt_year);
        if(flag_leap)
        {
            offset_date+=366.0;
        }
        else
        {
            offset_date+=365.0;
        }
        
    }
    
    if(isLeapYear(yint))
    {
        vec_accumulate_date[0]=0;
        vec_accumulate_date[1]=31;
        vec_accumulate_date[2]=60;
        vec_accumulate_date[3]=91;
        vec_accumulate_date[4]=121;
        vec_accumulate_date[5]=152;
        vec_accumulate_date[6]=182;
        vec_accumulate_date[7]=213;
        vec_accumulate_date[8]=244;
        vec_accumulate_date[9]=274;
        vec_accumulate_date[10]=305;
        vec_accumulate_date[11]=335;
    }
    else
    {
        vec_accumulate_date[0]=0;
        vec_accumulate_date[1]=31;
        vec_accumulate_date[2]=59;
        vec_accumulate_date[3]=90;
        vec_accumulate_date[4]=120;
        vec_accumulate_date[5]=151;
        vec_accumulate_date[6]=181;
        vec_accumulate_date[7]=212;
        vec_accumulate_date[8]=243;
        vec_accumulate_date[9]=273;
        vec_accumulate_date[10]=304;
        vec_accumulate_date[11]=334;
    }

    offset_date+=(double)vec_accumulate_date[mint-1] +(double)dint +Hint/24.0 +Mint/1440.0 +Sint/86400.0;
    
    return offset_date;
}

float get_dt(char *string_d0,char *string_d1)
{
    float out;
    double t0=get_datenum(string_d0);
    double t1=get_datenum(string_d1);
    double dt=t1-t0;

    out=(float)dt;

    return out;
}   //TODO: Implement it in more accurate and faster way

void getTimeStampStr(char *filename, char *str_timestamp)
{
    uint16_t index_last_slash;
    uint16_t cnt;
    //find the location of the last slash
    for(cnt=0;filename[cnt]!='\0';cnt++)
    {
        if(filename[cnt]=='/')
        {
            index_last_slash=cnt;
        }
    }

    for(cnt=0;cnt<14;cnt++)
    {
        str_timestamp[cnt]=filename[index_last_slash+1+cnt];
    }
    str_timestamp[14]='\0';


}

