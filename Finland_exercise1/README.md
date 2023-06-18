
The script **locate_one_event.sh** will locate the specified event relative to event H01 using the time-measurements in the file **../Finland_dataset/Finland_CC_times.txt** and the set of slowness vectors generated in the file **../Finland_dataset/Finland_ak135_slovecs.txt**.  
It assumes that the event is co-located with event H01 at the start and moves the event, with the slowness vectors unchanged, by gradient descent until the norm of the solution vector ceases to decrease.  

We type

```
./locate_one_event.sh     H02
```

to locate event H02 with respect to event H01.
