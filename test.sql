\timing
SET client_min_messages = 'DEBUG5';
SET log_min_messages = 'DEBUG5';
SET wal_level = 'minimal';

create extension if not exists cube;

begin transaction;
SELECT setseed(.43);

create table dataTable
as 
select x as group_id, y as entry_id, 1::decimal(19,6) as heavy0  from generate_series(1,1e2,1) x, generate_series(1,1e5,1) y;

create index idx on dataTable(group_id,entry_id) including (heavy0);

select pg_size_pretty(pg_relation_size('idx'));

rollback;
