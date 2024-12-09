#include<stdio.h>
#include<string.h>
#include<stdlib.h>
typedef struct scores{
        char id[15];
        char name[10];
        int mark[4];
        int sum;
        float ave;
        struct scores *next;
}score;

void sortdata(score **headp){
    score *head, *p, *last, *tail;
    head = (score*)malloc(sizeof(score));
    head->next = *headp;
    for(tail = NULL; head->next != tail; tail = p){
        for(last = head, p = head->next; p->next != tail; last = p, p = p->next){
            if(p->sum > p->next->sum){
                last->next = p->next;
                p->next = p->next->next;
                last->next->next = p;
                p = last->next;
            }
        }
    }
    *headp = head->next;
    free(head);
    
}

int main()
{
    score *p = NULL, *tail = NULL, *head = NULL, *move = NULL, *move1 = NULL, *tmove = NULL;
    
    int num, t, i;
    int sig, newscore;
    char check[15];
    int con;
    int flag = 0, flag2 = 0;
do{

    scanf("%d",&con);
    switch(con){
        case 1:
            scanf("%d",&num);
            t=0;
            while(t < num){
                p = (score *)malloc(sizeof(score));
                scanf("%s %s %d %d %d %d", p->id, p->name, &(p->mark[0]), &(p->mark[1]), &(p->mark[2]), &(p->mark[3]));
                p->sum = (p->mark[0] + p->mark[1] + p->mark[2] + p->mark[3]);
                if(head == NULL){
                    head = p;
                }else{
                    tail->next = p;
                }
                tail = p;
                t++;
            }
            sortdata(&head);
            
            
                
        case 2:    
            tmove=head;
            while(tmove!=NULL){
                printf("%s %s %d %d %d %d\n",tmove->id,tmove->name,tmove->mark[0],tmove->mark[1],tmove->mark[2],tmove->mark[3]);
                tmove=tmove->next;
                
            }
            break;
        case 3:
            tmove=head;
            scanf("%s %d %d",check,&sig,&newscore);
            while(tmove!=NULL){
                if(strcmp(check,tmove->id)==0){
                    tmove->mark[sig-1]=newscore;
                    break;
                }else{
                    tmove=tmove->next;
                }
                
            }
            break;
        case 4:
            tmove=head;
            while(tmove!=NULL){
            tmove->ave=1.0*(tmove->sum)/4;
            printf("%s %s %.2f\n",tmove->id,tmove->name,tmove->ave);
            tmove=tmove->next;
            }
            break;
        case 5:
            tmove=head;
            while(tmove!=NULL){
                printf("%s %s %d %.2f\n",tmove->id,tmove->name,tmove->sum,tmove->ave);
                tmove=tmove->next;
                
            }
            break;
            
    }

}while(con!=0);


return 0;
    
}