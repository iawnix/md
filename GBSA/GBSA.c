# include "stdio.h"
# include "stdlib.h"
# include "string.h"
# define SIGN1 "Total Energy Decomposition:"
# define SIGN2 "Sidechain Energy Decomposition:"
# define SIGN3 "Backbone Energy Decomposition:"

/*
Author zjh HENU
本质上是调用shell指令
*/

// 此函数用于找到对应的SIGN1
/*
 * part1 SIGN1 + 3 ~ SIGN2 - 1 
 * part2 SIGN2 + 3 ~ SIGN3 - 1
 * part3 SIGN3 + 3 ~ -1
 * */
int GetLineNum(char *filename,int i)
{
    FILE * fp = NULL;
    char Buf1[1024];
    char Buf2[1024];
    char CMD[1024];
    
    switch (i)
    {
        case 0 : sprintf(CMD,"grep -n '%s' %s | cut -f1 -d: ",SIGN1,filename); break;
        case 1 : sprintf(CMD,"grep -n '%s' %s | cut -f1 -d: ",SIGN2,filename); break;
        case 2 : sprintf(CMD,"grep -n '%s' %s | cut -f1 -d: ",SIGN3,filename); break;
        default : printf("error\n");
    }

    fp = popen(CMD,"r" );
    if (fp)
    {
        fgets(Buf1,1024,fp);
        sscanf(Buf1,"%s:%*s:",Buf2);

    }
    fclose(fp);
    
    return atoi(Buf2);
}


void Devide(int i, int j,char *filename,char *out)
{
    char CMD[1024];
	
	// 2023-04-13 zjh 修订残基号四位数导致与残基名相连，就是 变,
	if (j == -1)
    {
        sprintf(CMD," sed -n '%d,%sp' %s  | awk -F ',' '{printf (\"%%s,%%s,%%s,%%s,%%s,%%s \\n\",$1,$6,$9,$12,$15,$18)}' > %s ",i,"$",filename,out);
    }
    else
    {
        sprintf(CMD," sed -n '%d,%dp' %s | awk -F ',' '{printf (\"%%s,%%s,%%s,%%s,%%s,%%s \\n\",$1,$6,$9,$12,$15,$18)}' > %s ",i,j,filename,out); 
    }
    //printf("%s \n",CMD);
    if (system(CMD) != 0)
    {
        printf("ERROR \n");
    }
    
}

typedef struct  sortednum
{
    int numl[800];			// 储存符合要求的数据行号
    int len;                // 储存行号数目
}Sortn;

void stripN(char *str);

// 用于根据总能量进行排序，同时输出Tol.csv文件
void Sorted(Sortn *my,int num)
{
    
    FILE *fp = NULL;                           // 用于指向管道
    FILE *fp2 = NULL;                          // 用于指向输出的Tol.csv
    char Buf1[1024];
    char * strvar = NULL;                      // 用于指向分割字符串

    char name[16];                             // 储存Tol中的残基名
    char number[16];                              // 储存Tol中的残基号
    char E[64];                            // 储存Tol中的E

	//char *CMD = "cat ./Total.zjh |awk -F ',' '{printf (\"%s,%s\\n\",$1,$6)}' |sort -t ',' -g -k 2";
	char *CMD = "grep -n ',' ./Total.zjh | sort -t ',' -g -k 6 | cut -f1 -d:";
	char *CMD_ = "cat ./Total.zjh |awk -F ',' '{printf(\"%s,%s \\n\"),$1,$6}'";
    fp = popen(CMD,"r" );
    if (fp)
    {
        for (int i = 0; i < num; i++)
        {
            if (!feof(fp))
            {
                fgets(Buf1,1024,fp);
				//printf("%s",Buf1);
				my->numl[i] = atoi(Buf1);
				printf("%d \n",my->numl[i]);
                (my->len) += 1;
				//printf("{%d,%d,%p}\n",i,my->n,fp);
				//printf("{n = %d,foef(fp)=%d, i = %d, num = %d}\n",my->n,feof(fp),i,num);
            }
            else
            {
				break;
            }
        }

    }
    
	fclose(fp);

    // 处理Tol文件，主要是将name与num分开

	fp = popen(CMD_,"r");
    fp2 = fopen("Tol.csv","w");
    if (fp2 == NULL)
    {
        printf("Error: cannot create Tol.csv\n");
    }
 
    if (fp)
    {
        fgets(Buf1,1024,fp);
        while (!feof(fp))
        {
            // 读取name num部分
            strvar = strtok(Buf1,",");
            // 需要分割
            strncpy(name,strvar,3);
            name[3] = '\0';
            strcpy(number,strvar+3);
            // 读取后面的能量
            strvar = strtok(NULL,",");
            // 去除末尾的\n
            strcpy(E,strvar);
            stripN(E);
            // 写入Tol.csv文件
            fprintf(fp2,"%s , %s , %s \n",name,number,E);
            // 读取下一行
            fgets(Buf1,1024,fp);
        }
    }
}


typedef struct unit
{
    char name[16];
    char num[16];
    char vdW[64];
    char Ele[64];
    char PS[64];
    char noPS[64];
    char TOTAL[64];
    char Side[64];
    char Bone[64];
}UNIT;

// 此函数用于去除字符串末尾的换行符1
void stripN(char *str)
{
	char *find = strchr(str,'\n');
	if (find)
	{
		*find = '\0';
	}
}

// 此函数用于按照.分割读取的数据
// 需要注意删除末尾的\n
// sign 1 : Total 2 : Side 或者Bone
void Mydelim(char *str, int sign, UNIT *my)
{
	char * strvar;
	char var[1024];
	strvar = strtok(str,",");
    

	if (sign == 1)
	{
        // 需要分割开，前三个是name，剩下的是num
        strncpy(my->name,strvar,3);
        *(my->name + 3) = '\0';
        strcpy(my->num,strvar + 3);
		strvar = strtok(NULL, ",");
		sprintf(my->vdW,"%s",strvar); 
		strvar = strtok(NULL, ",");
		sprintf(my->Ele,"%s",strvar); 
		strvar = strtok(NULL, ",");
		sprintf(my->PS,"%s",strvar); 
		strvar = strtok(NULL, ",");
		sprintf(my->noPS,"%s",strvar); 
		strvar = strtok(NULL, ",");
		stripN(strvar);
		sprintf(my->TOTAL,"%s",strvar); 
	}
	
    else if(sign == 2)
	{
        
		while (strvar != NULL) 
		{
			sprintf(var,"%s",strvar); 
			strvar = strtok(NULL, ",");
			if (strvar == NULL)
			{
				stripN(var);
				sprintf(my->Side,"%s",var); 
                break;
			}
		}
	}
	else 
	{
		while (strvar != NULL) 
		{
			sprintf(var,"%s",strvar); 
			strvar = strtok(NULL, ",");
			if (strvar == NULL)
			{
				stripN(var);
				sprintf(my->Bone,"%s",var); 
			}
		}
    }
}

void outKEY(Sortn *my)
{
    UNIT myunit;
    char line1[] = "name,num,vdW,Ele,PS,nPS,TOTAL,Side,Bone";
    char CMD[1024];
    char Buf1[1024];
    FILE *fp =NULL;
    FILE *fp_ =NULL;

    fp_ = fopen("KEY.csv","w");

    if (fp_ == NULL)
    {
        printf("Error: cannot open KEY.csv\n");
    }
    fprintf(fp_,"%s\n",line1);

    for(int i = 0; i < my->len ;i ++)
    {
        sprintf(CMD,"sed -n '%dp' ./Total.zjh",my->numl[i]);
        fp = popen(CMD,"r" );
        //printf("i = %d, part = Total\n%s\n", i,CMD);
        if (fp)
        {
            fgets(Buf1,1024,fp);
			Mydelim(Buf1,1,&myunit);
        }
        fclose(fp);

        sprintf(CMD,"sed -n '%dp' ./Side.zjh",my->numl[i]);
        fp = popen(CMD,"r" );
        //printf("i = %d, part = Side\n", i);
        if (fp)
        {
            fgets(Buf1,1024,fp);
            //sscanf(Buf1,"%*s,%*s,%*s,%*s,%*s,%s \n",myunit.Side);
			Mydelim(Buf1,2,&myunit);
        }
        fclose(fp);

        sprintf(CMD,"sed -n '%dp' ./Bone.zjh",my->numl[i]);
        fp = popen(CMD,"r" );
        //printf("i = %d, part = Bone\n", i);
        if (fp)
        {
            fgets(Buf1,1024,fp);
            // sscanf(Buf1,"%*s,%*s,%*s,%*s,%*s,%s \n",myunit.Bone);
			Mydelim(Buf1,0,&myunit);
        }
        fclose(fp);
        fprintf(fp_,"%s , %s , %s , %s , %s , %s , %s , %s , %s \n",myunit.name,myunit.num,myunit.vdW,myunit.Ele,myunit.PS,myunit.noPS,myunit.TOTAL,myunit.Side,myunit.Bone);
    }
    fclose(fp_);
}

int main(int argc, char *argv[])
{
    char *filename = NULL;
    int num[3];
    int res;
    res = atoi(argv[2]);
    filename = argv[1];

    Sortn my;
    my.len = 0;
    for (int i = 0; i < 3; i++)
    {
        num[i] = GetLineNum(filename,i);
    }

    Devide(num[0]+3,num[1]-2,filename,"Total.zjh");
    Devide(num[1]+3,num[2]-2,filename,"Side.zjh");
    Devide(num[2]+3,-1,filename,"Bone.zjh");

    Sorted(&my,res);     
    
    outKEY(&my);
	system("rm Total.zjh Side.zjh Bone.zjh");

    return 0;
}

