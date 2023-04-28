from datetime import datetime
from datetime import timedelta

def get_current_date(state):
    init_date = str(*state["init_date"])
    seconds = int(*state["seconds"])
    dt = datetime.strptime(init_date, "%Y%m%d")
    current_date = dt + timedelta(seconds=seconds)

    state["current_date"] = int(current_date.strftime("%Y%m%d"))
    state["current_month"] = int(current_date.strftime("%m"))    
    state["current_day"] = int(current_date.strftime("%d"))
    state["current_year"] = int(current_date.strftime("%Y"))
