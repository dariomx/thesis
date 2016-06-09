select target_id,
       median(to_number(m.value)) as median,
            (max(to_number(m.value)) - min(to_number(m.value))) as range
from ems_bdasdev_metrics_val m
where m.target_id = :tid
  and m.metric_id = :mid
  and m.collection_time
      between to_date(:cstart, 'yyyy-mm-dd hh24:mi:ss')
          and to_date(:cend, 'yyyy-mm-dd hh24:mi:ss')
group by target_id
