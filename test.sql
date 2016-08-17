\timing
SET client_min_messages = 'DEBUG5';
SET log_min_messages = 'DEBUG5';
SET wal_level = 'minimal';

create extension if not exists cube;

begin transaction;
SELECT setseed(.43);
/*
create table dataTable_justindex
as 
select x as group_id, y as entry_id,
1::decimal(19,6) as heavy0, 2::decimal(19,6) as heavy1, 3::decimal(19,6) as heavy2, 4::decimal(19,6) as heavy3 
from generate_series(1,1e2,1) x, generate_series(1,1e6,1) y;

create index idx_justindex on dataTable_justindex(group_id,entry_id);

create table dataTable_oldinclude
as 
select x as group_id, y as entry_id,
1::decimal(19,6) as heavy0, 2::decimal(19,6) as heavy1, 3::decimal(19,6) as heavy2, 4::decimal(19,6) as heavy3 
from generate_series(1,1e2,1) x, generate_series(1,1e6,1) y;

create index idx_oldinclude on dataTable_oldinclude(group_id,entry_id,heavy0,heavy1,heavy2,heavy3) ;
	
*/

create table dataTable_including
as 
select x as group_id, y as entry_id,
1::decimal(19,6) as heavy0, 2::decimal(19,6) as heavy1, 3::decimal(19,6) as heavy2, 4::decimal(19,6) as heavy3
from generate_series(1,1e2,1) x, generate_series(1,1e6,1) y;

create index idx_including on dataTable_including(group_id,entry_id) including (heavy0,heavy1,heavy2,heavy3);


--testing adhoc pre-grouping arbitrary filter which cannot be indexed

--select group_id from dataTable_justindex where heavy0 + heavy1 + heavy2 + heavy3 > 10 group by group_id;
--select group_id from dataTable_oldinclude where heavy0 + heavy1 + heavy2 + heavy3 > 10 group by group_id;
select group_id from dataTable_including where heavy0 + heavy1 + heavy2 + heavy3 > 10 group by group_id;
	

--select pg_size_pretty(pg_relation_size('idx_justindex'));
--select pg_size_pretty(pg_relation_size('idx_oldinclude'));
select pg_size_pretty(pg_relation_size('idx_including'));


rollback;
