select metric_id
from   ems_bdasdev_metrics_def
where  metric_name = :metric_name
  and  column_name = :column_name
