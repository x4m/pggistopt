\timing

create extension if not exists cube;
begin transaction;
create table dataTable(c cube);
create index idx on dataTable using gist(c);

insert into dataTable(c) select cube(array[random(),random(),random()]) from generate_series(1,1e5,1);

create table queries(id int,l1 float,l2 float,l3 float, u1 float,u2 float, u3 float, q cube);
insert into queries(id,l1,l2,l3) select s,random(),random(),random() from generate_series(1,1e4,1) s;
update queries set q = cube(array[l1,l2,l3],array[l1+0.1,l2+0.1,l3+0.1]);

select id,(select count(*) from dataTable dt where dt.c<@q) from queries ;

rollback;
--drop index idx;
--drop table dataTable;
--drop table queries;
